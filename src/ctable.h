#ifndef CAMIL_CTABLE_H
#define CAMIL_CTABLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "logger.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


// The immutable structure every classification thread queries: a sorted array
// of 64 bit keys, a parallel array of values, and a directory that turns the
// top bits of a key straight into the short run of keys it could be in.
//
//   directory[2^dir_bits + 1]   offsets into keys/values
//   keys[nentries]              ascending
//   values[nentries]            species, level, flags
//
// Why sorted runs rather than an open addressed table. A hash table has to
// keep empty slots, and rounding its capacity to a power of two wastes up to
// as much again, which is how a few hundred million cores turned into four
// gigabytes. Here the arrays are exactly as long as there is data, and the
// directory costs half a byte per entry at the sizing used below. A core that
// several species share is simply a run of consecutive entries with the same
// key, which is also what lifts the old limit of three species per core.
//
// The cost is one extra dependent memory access: the directory has to be read
// before the keys can be. Classification hides it by prefetching both levels
// several cores ahead, which it can do because a read is parsed into all of
// its cores before any of them is looked up.
//
// Keys and values are separate arrays because most lookups miss. A miss reads
// one cache line of keys and never touches values at all.

// entries a directory bucket should span on average, one cache line of keys
#define CTABLE_BUCKET_TARGET 8

// bounds on the directory, so it neither degenerates nor dwarfs the data
#define CTABLE_MIN_DIR_BITS 8
#define CTABLE_MAX_DIR_BITS 26

// how many cores ahead the classification loop starts the directory fetch, and
// how many ahead it starts the key fetch that depends on it
#define CTABLE_PREFETCH_DISTANCE 12
#define CTABLE_PREFETCH_NEAR 6

// returned by ctable_find() when the key is absent
#define CTABLE_NOT_FOUND UINT64_MAX

struct ctable {
	uint32_t *directory;        // 2^dir_bits + 1 run offsets
	uint64_t *keys;             // nentries, ascending
	struct camil_value *values; // nentries, parallel to keys
	uint64_t nentries;
	uint32_t dir_bits;
};

// Allocates room for `nentries`, choosing the directory size. The caller fills
// keys and values in ascending key order and then calls ctable_index().
// Returns 0 on success, -1 on allocation failure.
int ctable_alloc(struct ctable *table, uint64_t nentries);

// Builds the directory from the keys, which must already be sorted.
void ctable_index(struct ctable *table);

// Releases the table.
void ctable_free(struct ctable *table);

// Approximate heap footprint in bytes.
uint64_t ctable_memory(const struct ctable *table);

// Number of entries stored.
static inline uint64_t ctable_size(const struct ctable *table) {
	return table->nentries;
}

// Reads and writes the three arrays. Both return 0 on success, -1 on an I/O
// error or a malformed table.
int ctable_write(const struct ctable *table, FILE *out);
int ctable_read(struct ctable *table, FILE *in);

// Checks a table that came from a file: keys ascending, directory consistent,
// species ids below `ngenomes`. Returns 0 when sound, -1 otherwise.
int ctable_validate(const struct ctable *table, uint32_t ngenomes);

// Bucket a key falls in.
static inline uint32_t ctable_bucket(const struct ctable *table, uint64_t key) {
	return (uint32_t)(key >> (64 - table->dir_bits));
}

// Starts the fetch of the directory entry for `bucket`. The first of the two
// stages a caller should run ahead of the lookup.
static inline void ctable_prefetch_dir(const struct ctable *table, uint32_t bucket) {
#if defined(__GNUC__)
	__builtin_prefetch(&table->directory[bucket], 0, 3);
#else
	(void)table;
	(void)bucket;
#endif
}

// Starts the fetch of the keys the bucket covers. Only useful once the
// directory entry itself has arrived, so callers run it one stage behind
// ctable_prefetch_dir().
static inline void ctable_prefetch_keys(const struct ctable *table, uint32_t bucket) {
#if defined(__GNUC__)
	__builtin_prefetch(&table->keys[table->directory[bucket]], 0, 3);
#else
	(void)table;
	(void)bucket;
#endif
}

// Finds `key`. Returns the index of its first entry and writes the run length
// to `count`, or CTABLE_NOT_FOUND when the key is absent. Lock free, and safe
// from any number of threads since the table is never written after building.
static inline uint64_t ctable_find(const struct ctable *table, uint64_t key, uint32_t *count) {
	uint32_t bucket = ctable_bucket(table, key);
	uint64_t i = table->directory[bucket];
	uint64_t end = table->directory[bucket + 1];
	uint64_t first;

	// The run is short by construction, so a scan beats a binary search and
	// the ascending order lets it stop early on a miss.
	while (i < end && table->keys[i] < key) {
		i++;
	}
	if (i == end || table->keys[i] != key) {
		return CTABLE_NOT_FOUND;
	}

	first = i;
	while (i < end && table->keys[i] == key) {
		i++;
	}
	*count = (uint32_t)(i - first);
	return first;
}

#ifdef __cplusplus
}
#endif

#endif
