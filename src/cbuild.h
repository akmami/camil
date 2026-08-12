#ifndef CAMIL_CBUILD_H
#define CAMIL_CBUILD_H


#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "ctable.h"
#include "logger.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdint.h>


// occupancy a shard is allowed to reach before it doubles
#define CBUILD_LOAD_FACTOR 0.70

// shards allocated per worker thread, and the range the result is clamped to
#define CBUILD_SHARDS_PER_THREAD 32u
#define CBUILD_MIN_SHARDS 32u
#define CBUILD_MAX_SHARDS 4096u

#define CBUILD_MIN_CAPACITY 16u

// one independently locked partition
struct cbuild_shard {
	struct camil_entry *entries;
	uint64_t capacity; // power of two
	uint64_t size;     // occupied slots
	uint64_t limit;    // grow when size reaches this
	pthread_mutex_t mutex;
	char padding[CAMIL_CACHE_LINE];
};

struct cbuild {
	struct cbuild_shard *shards;
	uint32_t nshards;    // power of two
	uint32_t shard_mask; // nshards - 1
};

// per thread scratch used by cbuild_insert_cores()
struct cbuild_scratch {
	lcp_label *labels; // labels grouped by shard
	size_t capacity;   // labels the buffer can hold
	uint32_t *counts;  // per shard label count
	uint32_t *offsets; // per shard write cursor
	uint32_t nshards;
};

// creates a builder with a shard count derived from threads and room for capacity_hint entries before the first rehash
int cbuild_init(struct cbuild *builder, int threads, uint64_t capacity_hint);

// releases every shard
void cbuild_free(struct cbuild *builder);

// prepares scratch space for one worker thread. Returns 0 on success
int cbuild_scratch_init(struct cbuild_scratch *scratch, const struct cbuild *builder);

// releases the scratch space
void cbuild_scratch_free(struct cbuild_scratch *scratch);

// records that every core in cores occurs in genome gid
// thread safe with respect to other threads calling it, provided each passes its own scratch
int cbuild_insert_cores(struct cbuild *builder, struct cbuild_scratch *scratch, const struct core *cores, int ncores, camil_gid gid);

// total number of entries currently held, overflowed ones included
uint64_t cbuild_size(const struct cbuild *builder);

// builds the immutable table: every entry that belongs to between 1 and max_share genomes is copied into table, the rest are dropped
// builder is left untouched and should be freed by the caller afterwards
int cbuild_freeze(const struct cbuild *builder, int max_share, struct ctable *table);

#ifdef __cplusplus
}
#endif

#endif
