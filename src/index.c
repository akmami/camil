#include "index.h"


// number of builder slots preallocated before the first rehash
#define CAMIL_INDEX_INITIAL_SLOTS (1ull << 20)

// counting semaphore over sequence bases, used to cap how much sequence data is queued for parsing at any moment
struct index_budget {
	pthread_mutex_t mutex;
	pthread_cond_t cv;
	uint64_t used;
	uint64_t cap;
};

// fixed set of scratch buffers handed out one per running job
struct scratch_pool {
	struct cbuild_scratch *items;
	int *available; // stack of free indices
	int nfree;
	int count;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
};

// shared state for one indexing run
struct index_builder {
	struct camil_index *index;
	struct cbuild table;
	struct tpool pool;
	struct index_budget budget;
	struct scratch_pool scratch;
	pthread_mutex_t stats_mutex; // guards the counters below
	uint64_t nseq;               // sequences parsed for the current genome
	uint64_t nbases;             // bases parsed for the current genome
	uint64_t ncores;             // cores produced for the current genome
	int failed;                  // set when a worker hits an unrecoverable error
};

// one queued sequence. The worker owns seq and frees it when done
struct index_job {
	struct index_builder *builder;
	char *seq;
	int len;
	camil_gid gid;
};

static void budget_init(struct index_budget *budget, uint64_t cap) {
	pthread_mutex_init(&budget->mutex, NULL);
	pthread_cond_init(&budget->cv, NULL);
	budget->used = 0;
	budget->cap = cap;
}

static void budget_destroy(struct index_budget *budget) {
	pthread_mutex_destroy(&budget->mutex);
	pthread_cond_destroy(&budget->cv);
}

// reserves amount bases, waiting while that would exceed the cap
static void budget_acquire(struct index_budget *budget, uint64_t amount) {
	pthread_mutex_lock(&budget->mutex);
	while (budget->used > 0 && budget->used + amount > budget->cap) {
		pthread_cond_wait(&budget->cv, &budget->mutex);
	}
	budget->used += amount;
	pthread_mutex_unlock(&budget->mutex);
}

static void budget_release(struct index_budget *budget, uint64_t amount) {
	pthread_mutex_lock(&budget->mutex);
	budget->used -= amount;
	pthread_cond_broadcast(&budget->cv);
	pthread_mutex_unlock(&budget->mutex);
}

static int scratch_pool_init(struct scratch_pool *pool, const struct cbuild *table, int count) {
	int i;

	memset(pool, 0, sizeof(*pool));
	if (count < 1) {
		count = 1;
	}

	pool->items = (struct cbuild_scratch *)calloc((size_t)count, sizeof(struct cbuild_scratch));
	pool->available = (int *)malloc((size_t)count * sizeof(int));
	if (pool->items == NULL || pool->available == NULL) {
		log_error("out of memory while allocating worker scratch space");
		free(pool->items);
		free(pool->available);
		memset(pool, 0, sizeof(*pool));
		return -1;
	}

	// initialized before the loop below can fail
	pthread_mutex_init(&pool->mutex, NULL);
	pthread_cond_init(&pool->cv, NULL);

	for (i = 0; i < count; i++) {
		if (cbuild_scratch_init(&pool->items[i], table) != 0) {
			pool->count = i;
			return -1;
		}
		pool->available[i] = i;
	}

	pool->count = count;
	pool->nfree = count;
	return 0;
}

static void scratch_pool_destroy(struct scratch_pool *pool) {
	int i;

	for (i = 0; i < pool->count; i++) {
		cbuild_scratch_free(&pool->items[i]);
	}
	// the lock exists whenever the arrays do, including after a partially failed setup
	if (pool->items != NULL) {
		pthread_mutex_destroy(&pool->mutex);
		pthread_cond_destroy(&pool->cv);
	}
	free(pool->items);
	free(pool->available);
	memset(pool, 0, sizeof(*pool));
}

static struct cbuild_scratch *scratch_acquire(struct scratch_pool *pool, int *slot) {
	pthread_mutex_lock(&pool->mutex);
	while (pool->nfree == 0) {
		pthread_cond_wait(&pool->cv, &pool->mutex);
	}
	*slot = pool->available[--pool->nfree];
	pthread_mutex_unlock(&pool->mutex);
	return &pool->items[*slot];
}

static void scratch_release(struct scratch_pool *pool, int slot) {
	pthread_mutex_lock(&pool->mutex);
	pool->available[pool->nfree++] = slot;
	pthread_cond_signal(&pool->cv);
	pthread_mutex_unlock(&pool->mutex);
}

// records a failure in a way the reading thread can observe.
static void index_fail(struct index_builder *builder) {
	pthread_mutex_lock(&builder->stats_mutex);
	builder->failed = 1;
	pthread_mutex_unlock(&builder->stats_mutex);
}

// reading thread polls this so it stops queueing work once the table can no longer accept it
static int index_failed(struct index_builder *builder) {
	int failed;

	pthread_mutex_lock(&builder->stats_mutex);
	failed = builder->failed;
	pthread_mutex_unlock(&builder->stats_mutex);
	return failed;
}

// parses one orientation of a sequence and folds every core label into the builder under gid
static uint64_t index_parse(struct index_builder *builder, struct cbuild_scratch *scratch, const char *seq, int len, camil_gid gid, int reverse_complement) {
	struct lps parsed;
	uint64_t ncores;

	if (reverse_complement) {
		init_lps2(&parsed, seq, len);
	} else {
		init_lps(&parsed, seq, len);
	}
	lps_deepen(&parsed, builder->index->lcp_level);

	if (cbuild_insert_cores(&builder->table, scratch, parsed.cores, parsed.size, gid) != 0) {
		index_fail(builder);
	}

	ncores = parsed.size > 0 ? (uint64_t)parsed.size : 0;
	free_lps(&parsed);
	return ncores;
}

// thread pool entry point: parse one queued sequence in both orientations
static void index_job_run(void *arg) {
	struct index_job *job = (struct index_job *)arg;
	struct index_builder *builder = job->builder;
	struct cbuild_scratch *scratch;
	uint64_t ncores;
	int slot;

	scratch = scratch_acquire(&builder->scratch, &slot);

	ncores = index_parse(builder, scratch, job->seq, job->len, job->gid, 0);
	// skip the second orientation once the table has stopped accepting insertions; there is nothing left for it to do
	if (builder->index->use_rc && !index_failed(builder)) {
		ncores += index_parse(builder, scratch, job->seq, job->len, job->gid, 1);
	}

	scratch_release(&builder->scratch, slot);

	free(job->seq);
	budget_release(&builder->budget, (uint64_t)job->len);

	pthread_mutex_lock(&builder->stats_mutex);
	builder->nseq++;
	builder->nbases += (uint64_t)job->len;
	builder->ncores += ncores;
	pthread_mutex_unlock(&builder->stats_mutex);

	free(job);
}

// reads one reference file and queues every sequence in it
static int index_add_genome(struct index_builder *builder, const char *path, camil_gid gid) {
	struct seqfile *file;
	struct seqrec rec;
	int ret;
	int status = 0;

	file = seq_open(path);
	if (file == NULL) {
		return -1;
	}

	pthread_mutex_lock(&builder->stats_mutex);
	builder->nseq = 0;
	builder->nbases = 0;
	builder->ncores = 0;
	pthread_mutex_unlock(&builder->stats_mutex);

	while ((ret = seq_read(file, &rec)) > 0) {
		struct index_job *job;
		char *copy;

		if (index_failed(builder)) {
			status = -1;
			break;
		}

		if (rec.len == 0) {
			continue;
		}

		if (rec.len > (size_t)INT_MAX) {
			log_warn("%s: sequence %s is %zu bases and exceeds the parser limit, skipped", path, rec.name, rec.len);
			continue;
		}

		copy = (char *)malloc(rec.len + 1);
		if (copy == NULL) {
			log_error("out of memory while reading %s", path);
			status = -1;
			break;
		}
		memcpy(copy, rec.seq, rec.len + 1);

		job = (struct index_job *)malloc(sizeof(struct index_job));
		if (job == NULL) {
			log_error("out of memory while reading %s", path);
			free(copy);
			status = -1;
			break;
		}
		job->builder = builder;
		job->seq = copy;
		job->len = (int)rec.len;
		job->gid = gid;

		// throttle the reader so that the queued sequences stay within budget
		budget_acquire(&builder->budget, rec.len);
		tpool_submit(&builder->pool, index_job_run, job);
	}

	if (ret < 0) {
		status = -1;
	}

	// finish the genome before moving on: the per genome statistics below assume one genome at a time
	tpool_wait(&builder->pool);
	seq_close(file);

	pthread_mutex_lock(&builder->stats_mutex);
	if (builder->failed) {
		status = -1;
	}
	pthread_mutex_unlock(&builder->stats_mutex);

	if (status == 0) {
		log_info("indexed %s (%s): %llu sequences, %llu bases, %llu cores", builder->index->names[gid], path, (unsigned long long)builder->nseq, (unsigned long long)builder->nbases, (unsigned long long)builder->ncores);
	}

	return status;
}

// index of the genome already reported under name, or -1 when it is free
static int name_owner(char *const *names, uint32_t count, const char *name) {
	uint32_t i;

	for (i = 0; i < count; i++) {
		if (names[i] != NULL && strcmp(names[i], name) == 0) {
			return (int)i;
		}
	}
	return -1;
}

// picks the label each genome is reported under. its short name when one was given, otherwise the file name with the directory and extension removed.
static int index_resolve_names(struct camil_index *index, const struct camil_genome *genomes) {
	uint32_t i;

	index->names = (char **)calloc(index->ngenomes, sizeof(char *));
	if (index->names == NULL) {
		log_error("out of memory while allocating genome names");
		return -1;
	}

	for (i = 0; i < index->ngenomes; i++) {
		char *candidate;
		int owner;

		candidate = genomes[i].name != NULL ? camil_strdup(genomes[i].name) : seq_basename(genomes[i].path);
		if (candidate == NULL) {
			log_error("out of memory while allocating genome names");
			return -1;
		}

		// requested name that is taken gives way to the file name
		owner = name_owner(index->names, i, candidate);
		if (owner >= 0 && genomes[i].name != NULL) {
			char *fallback = seq_basename(genomes[i].path);

			if (fallback == NULL) {
				log_error("out of memory while allocating genome names");
				free(candidate);
				return -1;
			}
			log_warn("the name '%s' is already used by %s, %s falls back to '%s'", candidate, genomes[owner].path, genomes[i].path, fallback);
			free(candidate);
			candidate = fallback;
		}

		// two files can still share a base name, so the last resort is a suffix
		owner = name_owner(index->names, i, candidate);
		if (owner >= 0) {
			char numbered[CAMIL_NAME_MAX];
			unsigned int suffix;

			for (suffix = 2; suffix <= index->ngenomes + 1; suffix++) {
				snprintf(numbered, sizeof(numbered), "%.*s.%u", (int)sizeof(numbered) - 12, candidate, suffix);
				if (name_owner(index->names, i, numbered) < 0) {
					break;
				}
			}
			log_warn("the name '%s' is already used by %s, %s is reported as '%s'", candidate, genomes[owner].path, genomes[i].path, numbered);
			free(candidate);
			candidate = camil_strdup(numbered);
			if (candidate == NULL) {
				log_error("out of memory while allocating genome names");
				return -1;
			}
		}

		index->names[i] = candidate;
		log_info("genome %u: %s (%s)", i, index->names[i], genomes[i].path);
	}

	return 0;
}

// releases everything an indexing run allocated apart from the index itself
static void index_builder_destroy(struct index_builder *builder) {
	tpool_destroy(&builder->pool);
	scratch_pool_destroy(&builder->scratch);
	cbuild_free(&builder->table);
	budget_destroy(&builder->budget);
	pthread_mutex_destroy(&builder->stats_mutex);
}

int camil_index_build(struct camil_index *index, const struct camil_genome *genomes, uint32_t ngenomes, int lcp_level, int max_share, int use_rc, int threads) {
	struct index_builder builder;
	uint32_t i;
	uint64_t collected;

	memset(index, 0, sizeof(*index));

	if (ngenomes == 0) {
		log_error("no reference genome was given");
		return -1;
	}
	if (ngenomes > CAMIL_MAX_GENOMES) {
		log_error("at most %u reference genomes can be indexed, %u were given", CAMIL_MAX_GENOMES, ngenomes);
		return -1;
	}
	if (lcp_level < 1) {
		log_error("the LCP level must be at least 1");
		return -1;
	}
	if (max_share < 1 || max_share > CAMIL_MAX_SHARE) {
		log_error("the sharing limit must be between 1 and %d", CAMIL_MAX_SHARE);
		return -1;
	}

	index->lcp_level = lcp_level;
	index->max_share = max_share;
	index->use_rc = use_rc != 0;
	index->ngenomes = ngenomes;

	if (index_resolve_names(index, genomes) != 0) {
		camil_index_destroy(index);
		return -1;
	}

	memset(&builder, 0, sizeof(builder));
	builder.index = index;
	pthread_mutex_init(&builder.stats_mutex, NULL);
	budget_init(&builder.budget, CAMIL_INDEX_INFLIGHT_BASES);

	if (cbuild_init(&builder.table, threads, CAMIL_INDEX_INITIAL_SLOTS) != 0 ||
	    scratch_pool_init(&builder.scratch, &builder.table, threads) != 0 ||
	    tpool_init(&builder.pool, threads, (size_t)(threads > 0 ? threads : 1) * 2) != 0) {
		index_builder_destroy(&builder);
		camil_index_destroy(index);
		return -1;
	}

	for (i = 0; i < ngenomes; i++) {
		if (index_add_genome(&builder, genomes[i].path, (camil_gid)i) != 0) {
			index_builder_destroy(&builder);
			camil_index_destroy(index);
			return -1;
		}
	}

	// every reference has been parsed; the mutable structure is no longer needed and is condensed into the layout classification wants
	collected = cbuild_size(&builder.table);
	if (cbuild_freeze(&builder.table, max_share, &index->table) != 0) {
		index_builder_destroy(&builder);
		camil_index_destroy(index);
		return -1;
	}
	index_builder_destroy(&builder);

	log_info("collected %llu distinct cores, %llu remain after dropping cores shared by more than %d genome(s)",
	         (unsigned long long)collected, (unsigned long long)ctable_size(&index->table), max_share);

	return 0;
}

int camil_index_save(const struct camil_index *index, const char *path) {
	FILE *out;
	uint32_t value;
	uint32_t i;

	out = fopen(path, "wb");
	if (out == NULL) {
		log_error("cannot write the index to %s", path);
		return -1;
	}

	if (fwrite(CAMIL_INDEX_MAGIC, 1, 8, out) != 8) {
		goto write_error;
	}

	value = CAMIL_INDEX_VERSION;
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}

	// recorded so that an index built with 64 bit labels is not silently read by a binary compiled with 32 bit labels
	value = (uint32_t)(sizeof(lcp_label) * 8);
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}

	value = (uint32_t)index->lcp_level;
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}
	value = (uint32_t)index->max_share;
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}
	value = (uint32_t)index->use_rc;
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}
	value = index->ngenomes;
	if (fwrite(&value, sizeof(value), 1, out) != 1) {
		goto write_error;
	}

	for (i = 0; i < index->ngenomes; i++) {
		uint32_t len = (uint32_t)strlen(index->names[i]);

		if (fwrite(&len, sizeof(len), 1, out) != 1 ||
		    fwrite(index->names[i], 1, len, out) != len) {
			goto write_error;
		}
	}

	if (ctable_write(&index->table, out) != 0) {
		goto write_error;
	}

	if (fclose(out) != 0) {
		log_error("cannot flush the index to %s", path);
		return -1;
	}

	log_info("wrote %llu cores to %s", (unsigned long long)ctable_size(&index->table), path);
	return 0;

write_error:
	log_error("failed while writing the index to %s", path);
	fclose(out);
	return -1;
}

int camil_index_load(struct camil_index *index, const char *path) {
	FILE *in;
	char magic[8];
	uint32_t version;
	uint32_t label_bits;
	uint32_t value;
	uint32_t g;

	memset(index, 0, sizeof(*index));

	in = fopen(path, "rb");
	if (in == NULL) {
		log_error("cannot open the index %s", path);
		return -1;
	}

	if (fread(magic, 1, 8, in) != 8 || memcmp(magic, CAMIL_INDEX_MAGIC, 8) != 0) {
		log_error("%s is not a camil index", path);
		goto read_error;
	}
	if (fread(&version, sizeof(version), 1, in) != 1 || version != CAMIL_INDEX_VERSION) {
		log_error("%s was written by an incompatible version of " CAMIL_NAME, path);
		goto read_error;
	}
	if (fread(&label_bits, sizeof(label_bits), 1, in) != 1 ||
	    label_bits != (uint32_t)(sizeof(lcp_label) * 8)) {
		log_error("%s stores %u bit core labels but this build uses %u bit labels", path, label_bits, (uint32_t)(sizeof(lcp_label) * 8));
		goto read_error;
	}

	if (fread(&value, sizeof(value), 1, in) != 1 || value < 1 || value > INT_MAX) {
		goto truncated;
	}
	index->lcp_level = (int)value;

	if (fread(&value, sizeof(value), 1, in) != 1 || value < 1 || value > CAMIL_MAX_SHARE) {
		goto truncated;
	}
	index->max_share = (int)value;

	if (fread(&value, sizeof(value), 1, in) != 1) {
		goto truncated;
	}
	index->use_rc = (int)value;

	if (fread(&value, sizeof(value), 1, in) != 1 || value == 0 || value > CAMIL_MAX_GENOMES) {
		goto truncated;
	}
	index->ngenomes = value;

	index->names = (char **)calloc(index->ngenomes, sizeof(char *));
	if (index->names == NULL) {
		log_error("out of memory while loading %s", path);
		goto read_error;
	}

	for (g = 0; g < index->ngenomes; g++) {
		uint32_t len;

		if (fread(&len, sizeof(len), 1, in) != 1 || len > (1u << 16)) {
			goto truncated;
		}
		index->names[g] = (char *)malloc(len + 1);
		if (index->names[g] == NULL) {
			log_error("out of memory while loading %s", path);
			goto read_error;
		}
		if (fread(index->names[g], 1, len, in) != len) {
			goto truncated;
		}
		index->names[g][len] = '\0';
	}

	if (ctable_read(&index->table, in) != 0) {
		goto truncated;
	}
	if (ctable_validate(&index->table, index->ngenomes) != 0) {
		log_error("%s is corrupted", path);
		goto read_error;
	}

	fclose(in);
	log_info("loaded %llu cores from %s", (unsigned long long)ctable_size(&index->table), path);
	return 0;

truncated:
	log_error("%s is truncated or corrupted", path);

read_error:
	fclose(in);
	camil_index_destroy(index);
	return -1;
}

void camil_index_destroy(struct camil_index *index) {
	uint32_t i;

	if (index->names != NULL) {
		for (i = 0; i < index->ngenomes; i++) {
			free(index->names[i]);
		}
		free(index->names);
	}
	ctable_free(&index->table);
	memset(index, 0, sizeof(*index));
}

const char *camil_index_genome_name(const struct camil_index *index, camil_gid gid) {
	if (index->names == NULL || gid >= index->ngenomes) {
		return "?";
	}
	return index->names[gid];
}

uint64_t camil_index_size(const struct camil_index *index) {
	return ctable_size(&index->table);
}

void camil_index_report(const struct camil_index *index) {
	uint32_t i;

	log_info("index: LCP level %d, sharing limit %d, reverse complement %s", index->lcp_level, index->max_share, index->use_rc ? "on" : "off");
	for (i = 0; i < index->ngenomes; i++) {
		log_info("  genome %u: %s", i, index->names[i]);
	}
	log_info("index: %llu cores, about %llu MiB in memory", (unsigned long long)ctable_size(&index->table), (unsigned long long)(ctable_memory(&index->table) / (1024ull * 1024ull)));

	camil_index_warn_collisions(index);
}

// warns when the index fills enough of the label space for chance matches to distort the result.
void camil_index_warn_collisions(const struct camil_index *index) {
	double space = ldexp(1.0, (int)(sizeof(lcp_label) * 8));
	double occupancy = (double)ctable_size(&index->table) / space;

	if (occupancy < CAMIL_COLLISION_WARN) {
		return;
	}

	log_warn("the index fills %.1f%% of the %u bit label space, so roughly %.1f%% of read cores will match by chance", 100.0 * occupancy, (unsigned)(sizeof(lcp_label) * 8), 100.0 * occupancy);
	log_warn("a read with %d cores then carries a spurious hit %.0f%% of the time; consider --min-ratio below 1.0 and --min-hits above 1", CAMIL_COLLISION_READ_CORES, 100.0 * (1.0 - pow(1.0 - occupancy, CAMIL_COLLISION_READ_CORES)));
}
