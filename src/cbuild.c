#include "cbuild.h"


static uint64_t cbuild_next_pow2(uint64_t value) {
	uint64_t result = CBUILD_MIN_CAPACITY;

	while (result < value) {
		result <<= 1;
	}
	return result;
}

// shard a label belongs to
static inline uint32_t cbuild_shard_of(const struct cbuild *builder, uint64_t hash) {
	return (uint32_t)((hash >> 40) & builder->shard_mask);
}

// reinserts an entry during a rehash
static void shard_place(struct camil_entry *entries, uint64_t mask, const struct camil_entry *entry) {
	uint64_t idx = camil_hash(entry->label) & mask;

	while (entries[idx].ngenomes != 0) {
		idx = (idx + 1) & mask;
	}
	entries[idx] = *entry;
}

// doubles a shard
static int shard_grow(struct cbuild_shard *shard) {
	uint64_t capacity = shard->capacity * 2;
	struct camil_entry *entries;
	uint64_t i;

	entries = (struct camil_entry *)calloc(capacity, sizeof(struct camil_entry));
	if (entries == NULL) {
		log_error("index builder: out of memory while growing a shard to %llu slots", (unsigned long long)capacity);
		return -1;
	}

	for (i = 0; i < shard->capacity; i++) {
		if (shard->entries[i].ngenomes != 0) {
			shard_place(entries, capacity - 1, &shard->entries[i]);
		}
	}

	free(shard->entries);
	shard->entries = entries;
	shard->capacity = capacity;
	shard->limit = (uint64_t)((double)capacity * CBUILD_LOAD_FACTOR);
	return 0;
}

// inserts one label into a shard the caller has already locked
static int shard_insert(struct cbuild_shard *shard, lcp_label label, camil_gid gid) {
	uint64_t mask = shard->capacity - 1;
	uint64_t idx = camil_hash(label) & mask;
	struct camil_entry *entry;
	uint8_t i;

	if (shard->size + 1 >= shard->capacity) {
		log_error("index builder: a shard is full and could not be grown");
		return -1;
	}

	while (shard->entries[idx].ngenomes != 0) {
		if (shard->entries[idx].label == label) {
			entry = &shard->entries[idx];

			// already recorded for this genome: a repeat within the genome
			for (i = 0; i < entry->ngenomes && i < CAMIL_MAX_SHARE; i++) {
				if (entry->gids[i] == gid) {
					return 0;
				}
			}
			if (entry->ngenomes < CAMIL_MAX_SHARE) {
				entry->gids[entry->ngenomes++] = gid;
			} else {
				// present in more genomes than an entry can hold
				// it cannot discriminate and is marked for removal at freeze time
				entry->ngenomes = CAMIL_OVERFLOW;
			}
			return 0;
		}
		idx = (idx + 1) & mask;
	}

	entry = &shard->entries[idx];
	entry->label = label;
	entry->gids[0] = gid;
	entry->ngenomes = 1;
	shard->size++;

	if (shard->size >= shard->limit) {
		return shard_grow(shard);
	}
	return 0;
}

int cbuild_init(struct cbuild *builder, int threads, uint64_t capacity_hint) {
	uint64_t nshards = (uint64_t)(threads > 0 ? threads : 1) * CBUILD_SHARDS_PER_THREAD;
	uint64_t per_shard;
	uint32_t i;

	memset(builder, 0, sizeof(*builder));

	if (nshards < CBUILD_MIN_SHARDS) nshards = CBUILD_MIN_SHARDS;
	if (nshards > CBUILD_MAX_SHARDS) nshards = CBUILD_MAX_SHARDS;
	nshards = cbuild_next_pow2(nshards);

	builder->shards = (struct cbuild_shard *)calloc(nshards, sizeof(struct cbuild_shard));
	if (builder->shards == NULL) {
		log_error("index builder: out of memory while allocating %llu shards", (unsigned long long)nshards);
		return -1;
	}
	builder->nshards = (uint32_t)nshards;
	builder->shard_mask = (uint32_t)(nshards - 1);

	per_shard = cbuild_next_pow2((uint64_t)((double)(capacity_hint / nshards + 1) / CBUILD_LOAD_FACTOR) + 1);

	for (i = 0; i < builder->nshards; i++) {
		struct cbuild_shard *shard = &builder->shards[i];

		shard->entries = (struct camil_entry *)calloc(per_shard, sizeof(struct camil_entry));
		if (shard->entries == NULL) {
			log_error("index builder: out of memory while allocating shard %u", i);
			cbuild_free(builder);
			return -1;
		}
		shard->capacity = per_shard;
		shard->limit = (uint64_t)((double)per_shard * CBUILD_LOAD_FACTOR);
		pthread_mutex_init(&shard->mutex, NULL);
	}

	return 0;
}

void cbuild_free(struct cbuild *builder) {
	uint32_t i;

	if (builder->shards == NULL) {
		return;
	}
	for (i = 0; i < builder->nshards; i++) {
		free(builder->shards[i].entries);
		if (builder->shards[i].capacity != 0) {
			pthread_mutex_destroy(&builder->shards[i].mutex);
		}
	}
	free(builder->shards);
	memset(builder, 0, sizeof(*builder));
}

int cbuild_scratch_init(struct cbuild_scratch *scratch, const struct cbuild *builder) {
	memset(scratch, 0, sizeof(*scratch));

	scratch->counts = (uint32_t *)calloc(builder->nshards, sizeof(uint32_t));
	scratch->offsets = (uint32_t *)calloc(builder->nshards, sizeof(uint32_t));
	if (scratch->counts == NULL || scratch->offsets == NULL) {
		log_error("index builder: out of memory while allocating worker scratch space");
		cbuild_scratch_free(scratch);
		return -1;
	}
	scratch->nshards = builder->nshards;
	return 0;
}

void cbuild_scratch_free(struct cbuild_scratch *scratch) {
	free(scratch->labels);
	free(scratch->counts);
	free(scratch->offsets);
	memset(scratch, 0, sizeof(*scratch));
}

// makes sure the scratch label buffer holds at least needed labels
static int scratch_reserve(struct cbuild_scratch *scratch, size_t needed) {
	lcp_label *grown;

	if (needed <= scratch->capacity) {
		return 0;
	}
	grown = (lcp_label *)realloc(scratch->labels, needed * sizeof(lcp_label));
	if (grown == NULL) {
		log_error("index builder: out of memory while grouping %zu cores", needed);
		return -1;
	}
	scratch->labels = grown;
	scratch->capacity = needed;
	return 0;
}

int cbuild_insert_cores(struct cbuild *builder, struct cbuild_scratch *scratch, const struct core *cores, int ncores, camil_gid gid) {
	uint32_t shard_id;
	uint32_t cursor = 0;
	int i;

	if (ncores <= 0) return 0;
	if (scratch_reserve(scratch, (size_t)ncores) != 0) return -1;

	memset(scratch->counts, 0, builder->nshards * sizeof(uint32_t));

	// pass one: how many labels land in each shard
	for (i = 0; i < ncores; i++) {
		scratch->counts[cbuild_shard_of(builder, camil_hash(cores[i].label))]++;
	}

	// turn the counts into write cursors, one contiguous run per shard
	for (shard_id = 0; shard_id < builder->nshards; shard_id++) {
		scratch->offsets[shard_id] = cursor;
		cursor += scratch->counts[shard_id];
	}

	// pass two: scatter the labels into their runs, thread local, no locks
	for (i = 0; i < ncores; i++) {
		lcp_label label = cores[i].label;

		scratch->labels[scratch->offsets[cbuild_shard_of(builder, camil_hash(label))]++] = label;
	}

	// pass three: one lock per non empty shard, then a run of insertions that all hit the same shard and stay in cache
	cursor = 0;
	for (shard_id = 0; shard_id < builder->nshards; shard_id++) {
		struct cbuild_shard *shard = &builder->shards[shard_id];
		uint32_t count = scratch->counts[shard_id];
		uint32_t k;
		int status = 0;

		if (count == 0) {
			continue;
		}

		pthread_mutex_lock(&shard->mutex);
		for (k = 0; k < count; k++) {
			if (shard_insert(shard, scratch->labels[cursor + k], gid) != 0) {
				status = -1;
				break;
			}
		}
		pthread_mutex_unlock(&shard->mutex);

		if (status != 0) {
			return -1;
		}
		cursor += count;
	}

	return 0;
}

uint64_t cbuild_size(const struct cbuild *builder) {
	uint64_t total = 0;
	uint32_t i;

	for (i = 0; i < builder->nshards; i++) {
		total += builder->shards[i].size;
	}
	return total;
}

int cbuild_freeze(const struct cbuild *builder, int max_share, struct ctable *table) {
	uint64_t kept = 0;
	uint32_t s;
	uint64_t i;

	// count the survivors first so that the frozen table is allocated once at exactly the right size
	for (s = 0; s < builder->nshards; s++) {
		const struct cbuild_shard *shard = &builder->shards[s];

		for (i = 0; i < shard->capacity; i++) {
			uint8_t n = shard->entries[i].ngenomes;

			if (n > 0 && n <= max_share) {
				kept++;
			}
		}
	}

	if (ctable_init(table, kept) != 0) {
		return -1;
	}

	// Shards are visited in order, so the frozen table is filled in a
	// deterministic order and two runs over the same references produce
	// byte identical index files
	for (s = 0; s < builder->nshards; s++) {
		const struct cbuild_shard *shard = &builder->shards[s];

		for (i = 0; i < shard->capacity; i++) {
			uint8_t n = shard->entries[i].ngenomes;

			if (n > 0 && n <= max_share && ctable_insert(table, &shard->entries[i]) != 0) {
				ctable_free(table);
				return -1;
			}
		}
	}

	return 0;
}
