#ifndef CAMIL_CBUILD_H
#define CAMIL_CBUILD_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "ctable.h"
#include "logger.h"
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


// Index construction, arranged so that no two threads ever touch the same
// memory until the very last prefix sum.
//
// Collection. A worker parsing a sequence appends one record per core to its
// own sink. Each sink keeps a separate chunk list per bucket, where the bucket
// is the top bits of the key, so appending is a pointer bump with no lock, no
// atomic and no shared cache line.
//
// Reduction. A bucket is then owned by exactly one thread, which gathers that
// bucket's chunks from every sink, sorts them, collapses repeats and applies
// the sharing limit. Buckets are disjoint key ranges, so the sorted buckets
// concatenated in order are globally sorted: no merge step and no heap.
//
// Emission. Surviving counts are prefix summed, which is the one point where
// threads meet, and each bucket copies itself into the final arrays.
//
// The sort is a radix sort. A comparison sort over a few hundred million
// records through a function pointer costs minutes, and the keys are fixed
// width, so there is nothing to gain by comparing.

// bytes per chunk. Small enough that the slack across buckets and threads
// stays in the tens of megabytes, large enough that allocation is rare.
#define CBUILD_BLOCK_BYTES 32768u

// bounds on the bucket count, which is also the unit of parallelism in the
// reduction, so it wants to be several times the thread count
#define CBUILD_MIN_BUCKETS 256u
#define CBUILD_MAX_BUCKETS 4096u
#define CBUILD_BUCKETS_PER_THREAD 8u

// records below this go through insertion sort instead of a radix pass
#define CBUILD_SMALL_SORT 48u

#if defined(__GNUC__) || defined(__clang__)
#define CAMIL_PACKED __attribute__((__packed__))
#else
#define CAMIL_PACKED
#endif

// One core occurrence. Packed to twelve bytes rather than padded to sixteen,
// which is a quarter of the peak build memory; x86 and ARM both load the
// unaligned key without penalty worth measuring.
struct CAMIL_PACKED cbuild_rec {
	uint64_t key;
	camil_sid sid;
	uint8_t level;
	uint8_t flags;
};

struct cbuild_block {
	struct cbuild_block *next;
	uint32_t n;
	struct cbuild_rec recs[1]; // CBUILD_BLOCK_RECS of them, sized at allocation
};

// a bucket's chunk list inside one sink
struct cbuild_bucket {
	struct cbuild_block *head;
	struct cbuild_block *tail;
	uint64_t count;
};

// per thread collection point, handed out one at a time to running jobs
struct cbuild_sink {
	struct cbuild_bucket *buckets;
};

struct cbuild {
	struct cbuild_sink *sinks;
	int nsinks;
	uint32_t nbuckets;    // power of two
	uint32_t bucket_bits;
	int failed;           // set when a sink could not allocate
	pthread_mutex_t lock; // guards `failed` only
};

// Creates the sinks, one per worker. Returns 0 on success, -1 otherwise.
int cbuild_init(struct cbuild *builder, int threads);

// Releases every chunk still held.
void cbuild_free(struct cbuild *builder);

// Appends one core occurrence to a sink. Thread safe as long as each thread
// passes a sink of its own. Returns 0 on success, -1 when out of memory.
int cbuild_push(struct cbuild *builder, struct cbuild_sink *sink, uint64_t key, camil_sid sid, uint8_t level);

// Total occurrences collected so far, repeats included.
uint64_t cbuild_count(const struct cbuild *builder);

// Sorts, collapses repeats, drops cores shared by more than `max_share`
// species, and builds the final table. `threads` bucket reducers run in
// parallel. Returns 0 on success, -1 on failure.
int cbuild_freeze(struct cbuild *builder, uint32_t max_share, int threads, struct ctable *table);

#ifdef __cplusplus
}
#endif

#endif
