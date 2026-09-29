#include "cbuild.h"
#include "tpool.h"
#include <assert.h> // ksort.h asserts its digit width without including this
#include "ksort.h"  // from deps


// Sorting the collected records by key, using klib's radix sort.
//
// It is a most significant digit sort that permutes in place, so unlike a
// least significant digit sort it needs no second array to bounce between,
// which halves what a bucket costs while it is being reduced. Records below
// RS_MIN_SIZE fall through to an insertion sort inside the macro.
//
// The top bucket_bits of every key in a bucket are equal, so the first pass or
// two find a single non empty digit and simply recurse; that costs one linear
// scan each and is not worth special casing.
#define cbuild_rec_key(rec) ((rec).key)
KRADIX_SORT_INIT(cbuild, struct cbuild_rec, cbuild_rec_key, sizeof(uint64_t))


// records that fit in one chunk
#define CBUILD_BLOCK_RECS ((CBUILD_BLOCK_BYTES - sizeof(struct cbuild_block *) - sizeof(uint32_t)) / sizeof(struct cbuild_rec))

// what one bucket reduces to, filled by its owning thread and read by nobody else
struct cbuild_out {
	struct cbuild_rec *recs;
	uint64_t n;
};

// state shared by the reduction jobs. Only `failed` is written by more than
// one thread, and only under the builder lock.
struct cbuild_reduce {
	struct cbuild *builder;
	struct cbuild_out *out;
	uint32_t max_share;
	int failed;
};

struct cbuild_job {
	struct cbuild_reduce *reduce;
	uint32_t bucket;
};

static void cbuild_fail(struct cbuild *builder) {
	pthread_mutex_lock(&builder->lock);
	builder->failed = 1;
	pthread_mutex_unlock(&builder->lock);
}

static uint32_t cbuild_pow2(uint32_t value) {
	uint32_t result = 1;

	while (result < value) {
		result <<= 1;
	}
	return result;
}

int cbuild_init(struct cbuild *builder, int threads) {
	uint32_t nbuckets;
	int i;

	memset(builder, 0, sizeof(*builder));

	if (threads < 1) {
		threads = 1;
	}

	nbuckets = cbuild_pow2((uint32_t)threads * CBUILD_BUCKETS_PER_THREAD);
	if (nbuckets < CBUILD_MIN_BUCKETS) {
		nbuckets = CBUILD_MIN_BUCKETS;
	}
	if (nbuckets > CBUILD_MAX_BUCKETS) {
		nbuckets = CBUILD_MAX_BUCKETS;
	}

	builder->nbuckets = nbuckets;
	builder->bucket_bits = 0;
	while (((uint32_t)1 << builder->bucket_bits) < nbuckets) {
		builder->bucket_bits++;
	}

	builder->nsinks = threads;
	builder->sinks = (struct cbuild_sink *)calloc((size_t)threads, sizeof(struct cbuild_sink));
	if (builder->sinks == NULL) {
		log_error("index builder: out of memory while allocating sinks");
		return -1;
	}

	for (i = 0; i < threads; i++) {
		builder->sinks[i].buckets = (struct cbuild_bucket *)calloc(nbuckets, sizeof(struct cbuild_bucket));
		if (builder->sinks[i].buckets == NULL) {
			log_error("index builder: out of memory while allocating sink %d", i);
			cbuild_free(builder);
			return -1;
		}
	}

	pthread_mutex_init(&builder->lock, NULL);
	return 0;
}

void cbuild_free(struct cbuild *builder) {
	int i;
	uint32_t b;

	if (builder->sinks != NULL) {
		for (i = 0; i < builder->nsinks; i++) {
			if (builder->sinks[i].buckets == NULL) {
				continue;
			}
			for (b = 0; b < builder->nbuckets; b++) {
				struct cbuild_block *block = builder->sinks[i].buckets[b].head;

				while (block != NULL) {
					struct cbuild_block *next = block->next;

					free(block);
					block = next;
				}
			}
			free(builder->sinks[i].buckets);
		}
		free(builder->sinks);
		pthread_mutex_destroy(&builder->lock);
	}
	memset(builder, 0, sizeof(*builder));
}

int cbuild_push(struct cbuild *builder, struct cbuild_sink *sink, uint64_t key, camil_sid sid) {
	uint32_t b = (uint32_t)(key >> (64 - builder->bucket_bits));
	struct cbuild_bucket *bucket = &sink->buckets[b];
	struct cbuild_rec *rec;

	if (bucket->tail == NULL || bucket->tail->n == CBUILD_BLOCK_RECS) {
		struct cbuild_block *block = (struct cbuild_block *)malloc(sizeof(struct cbuild_block) + (CBUILD_BLOCK_RECS - 1) * sizeof(struct cbuild_rec));

		if (block == NULL) {
			log_error("index builder: out of memory while collecting cores");
			cbuild_fail(builder);
			return -1;
		}
		block->next = NULL;
		block->n = 0;
		if (bucket->tail == NULL) {
			bucket->head = block;
		} else {
			bucket->tail->next = block;
		}
		bucket->tail = block;
	}

	rec = &bucket->tail->recs[bucket->tail->n++];
	rec->key = key;
	rec->sid = sid;
	bucket->count++;
	return 0;
}

uint64_t cbuild_count(const struct cbuild *builder) {
	uint64_t total = 0;
	int i;
	uint32_t b;

	for (i = 0; i < builder->nsinks; i++) {
		for (b = 0; b < builder->nbuckets; b++) {
			total += builder->sinks[i].buckets[b].count;
		}
	}
	return total;
}

// Sorts the species ids of one run of equal keys, so that repeats of the same
// species end up adjacent. Runs are at most one entry per species and usually
// one entry long.
static void cbuild_sort_run(struct cbuild_rec *recs, uint64_t n) {
	uint64_t i;

	for (i = 1; i < n; i++) {
		struct cbuild_rec value = recs[i];
		uint64_t j = i;

		while (j > 0 && recs[j - 1].sid > value.sid) {
			recs[j] = recs[j - 1];
			j--;
		}
		recs[j] = value;
	}
}

// Collapses a sorted bucket in place: repeats of the same core in the same
// species become one record, and a core in more species than the limit allows
// is dropped whole. Returns the number of records kept.
static uint64_t cbuild_collapse(struct cbuild_rec *recs, uint64_t n, uint32_t max_share) {
	uint64_t out = 0;
	uint64_t i = 0;

	while (i < n) {
		uint64_t run_end = i;
		uint64_t distinct = 0;
		uint64_t j;

		while (run_end < n && recs[run_end].key == recs[i].key) {
			run_end++;
		}

		cbuild_sort_run(&recs[i], run_end - i);

		// Count the distinct species in the run, compacting the duplicates to
		// the front of it as we go.
		for (j = i; j < run_end; j++) {
			if (j == i || recs[j].sid != recs[j - 1].sid) {
				recs[i + distinct] = recs[j];
				distinct++;
			}
		}

		if (distinct <= max_share) {
			// memmove, not memcpy: the survivors may overlap their destination
			// when nothing earlier was dropped.
			memmove(&recs[out], &recs[i], distinct * sizeof(struct cbuild_rec));
			out += distinct;
		}

		i = run_end;
	}

	return out;
}

// Reduces one bucket: gather from every sink, sort, collapse, keep the result.
static void cbuild_reduce_bucket(void *arg) {
	struct cbuild_job *job = (struct cbuild_job *)arg;
	struct cbuild_reduce *reduce = job->reduce;
	struct cbuild *builder = reduce->builder;
	struct cbuild_rec *recs;
	uint64_t total = 0;
	uint64_t filled = 0;
	int i;

	for (i = 0; i < builder->nsinks; i++) {
		total += builder->sinks[i].buckets[job->bucket].count;
	}
	if (total == 0) {
		return;
	}

	recs = (struct cbuild_rec *)malloc(total * sizeof(struct cbuild_rec));
	if (recs == NULL) {
		log_error("index builder: out of memory while reducing a bucket of %llu cores", (unsigned long long)total);
		cbuild_fail(builder);
		return;
	}

	// Chunks are released as they are consumed, so the peak is the collected
	// data rather than the collected data plus a copy of it.
	for (i = 0; i < builder->nsinks; i++) {
		struct cbuild_bucket *bucket = &builder->sinks[i].buckets[job->bucket];
		struct cbuild_block *block = bucket->head;

		while (block != NULL) {
			struct cbuild_block *next = block->next;

			memcpy(&recs[filled], block->recs, block->n * sizeof(struct cbuild_rec));
			filled += block->n;
			free(block);
			block = next;
		}
		bucket->head = NULL;
		bucket->tail = NULL;
		bucket->count = 0;
	}

	radix_sort_cbuild(recs, recs + total);

	reduce->out[job->bucket].n = cbuild_collapse(recs, total, reduce->max_share);
	reduce->out[job->bucket].recs = recs;
}

// Runs the parallel reduction over every bucket. On success *out holds one
// reduced, sorted record array per bucket (NULL and 0 for an empty bucket) that
// the caller must free; the collected chunks are consumed either way.
static int cbuild_reduce(struct cbuild *builder, uint32_t max_share, int threads, struct cbuild_out **out) {
	struct cbuild_reduce reduce;
	struct cbuild_job *jobs;
	struct tpool pool;
	uint32_t b;
	int status;

	*out = NULL;

	reduce.builder = builder;
	reduce.max_share = max_share;
	reduce.failed = 0;
	reduce.out = (struct cbuild_out *)calloc(builder->nbuckets, sizeof(struct cbuild_out));
	jobs = (struct cbuild_job *)calloc(builder->nbuckets, sizeof(struct cbuild_job));
	if (reduce.out == NULL || jobs == NULL) {
		log_error("index builder: out of memory while starting the reduction");
		free(reduce.out);
		free(jobs);
		return -1;
	}

	if (tpool_init(&pool, threads, (size_t)builder->nbuckets) != 0) {
		free(reduce.out);
		free(jobs);
		return -1;
	}

	for (b = 0; b < builder->nbuckets; b++) {
		jobs[b].reduce = &reduce;
		jobs[b].bucket = b;
		tpool_submit(&pool, cbuild_reduce_bucket, &jobs[b]);
	}
	tpool_wait(&pool);
	tpool_destroy(&pool);
	free(jobs);

	pthread_mutex_lock(&builder->lock);
	status = builder->failed ? -1 : 0;
	pthread_mutex_unlock(&builder->lock);

	if (status != 0) {
		for (b = 0; b < builder->nbuckets; b++) {
			free(reduce.out[b].recs);
		}
		free(reduce.out);
		return -1;
	}

	*out = reduce.out;
	return 0;
}

static void cbuild_out_free(struct cbuild_out *out, uint32_t nbuckets) {
	uint32_t b;

	if (out == NULL) {
		return;
	}
	for (b = 0; b < nbuckets; b++) {
		free(out[b].recs);
	}
	free(out);
}

int cbuild_freeze(struct cbuild *builder, uint32_t max_share, int threads, struct ctable *table) {
	struct cbuild_out *out;
	uint64_t nentries = 0;
	uint64_t written = 0;
	uint32_t b;

	memset(table, 0, sizeof(*table));

	if (cbuild_reduce(builder, max_share, threads, &out) != 0) {
		return -1;
	}

	for (b = 0; b < builder->nbuckets; b++) {
		nentries += out[b].n;
	}

	if (ctable_alloc(table, nentries) != 0) {
		cbuild_out_free(out, builder->nbuckets);
		return -1;
	}

	// Buckets are disjoint key ranges reduced in order, so copying them one
	// after another leaves the whole array sorted. Each bucket is released as
	// soon as it has been copied, so the peak is well below twice the data.
	for (b = 0; b < builder->nbuckets; b++) {
		uint64_t i;

		for (i = 0; i < out[b].n; i++) {
			table->keys[written] = out[b].recs[i].key;
			table->values[written].sid = out[b].recs[i].sid;
			written++;
		}
		free(out[b].recs);
		out[b].recs = NULL;
	}
	free(out);

	ctable_index(table);
	return 0;
}

int cbuild_freeze_full(struct cbuild *builder, int threads, struct cfull *full) {
	struct cbuild_out *out;
	uint64_t ncores = 0;
	uint64_t npairs = 0;
	uint64_t cursor = 0;
	uint32_t b;

	if (cbuild_reduce(builder, CBUILD_KEEP_ALL, threads, &out) != 0) {
		return -1;
	}

	// One pass to size the map and the species array: a run of equal keys is
	// one core, and every record in it is one (core, species) pair.
	for (b = 0; b < builder->nbuckets; b++) {
		uint64_t i;

		for (i = 0; i < out[b].n; i++) {
			if (i == 0 || out[b].recs[i].key != out[b].recs[i - 1].key) {
				ncores++;
			}
		}
		npairs += out[b].n;
	}

	if (cfull_alloc(full, ncores, npairs) != 0) {
		cbuild_out_free(out, builder->nbuckets);
		return -1;
	}

	// Every core's species ids go to the flat array in one contiguous run,
	// ascending, and the map points at it. Species within a run come out of
	// the collapse already sorted.
	for (b = 0; b < builder->nbuckets; b++) {
		uint64_t i = 0;

		while (i < out[b].n) {
			uint64_t key = out[b].recs[i].key;
			uint64_t begin = cursor;
			uint64_t j = i;

			while (j < out[b].n && out[b].recs[j].key == key) {
				full->sids[cursor++] = out[b].recs[j].sid;
				j++;
			}
			if (cfull_put(full, key, begin, (uint32_t)(j - i)) != 0) {
				cbuild_out_free(out, builder->nbuckets);
				return -1;
			}
			i = j;
		}
		free(out[b].recs);
		out[b].recs = NULL;
	}
	free(out);

	full->nsids = cursor;
	return 0;
}
