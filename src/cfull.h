// ---------------------------------------------------------------------------
// cfull.h -- the full core index: every core of every species, kept so that
// any subset of the species can be carved out later without reparsing.
//
// Layout
// ------
// A hash table (khashl) from the 64 bit core key to a packed value, and one
// flat array of species ids:
//
//   map:   key -> (offset, count)          one entry per distinct core
//   sids:  species ids, ascending per run  one entry per (core, species) pair
//
// The species a core occurs in are the run sids[offset .. offset+count). The
// value packs the offset into 40 bits and the count into 24, which is what
// bounds an index at a trillion pairs and sixteen million species.
//
// Nothing is filtered when the full index is built: a core seen in every
// genome is stored with every genome. The sharing limit is applied only when
// a subset is extracted, and it counts species inside the subset, so a core
// that is unique among the chosen species is kept even if the whole index
// shares it widely.
//
// The table is an ensemble of khashl sub tables, each owning the keys whose
// hash ends in its number. That keeps every sub table well inside the 32 bit
// addressing of one khashl table while the whole may hold billions of cores,
// and it makes the on disk form a plain dump of the sub tables' arrays, so
// loading is a read with no rehash.
//
// On disk format
// --------------
// Host byte order throughout.
//
//   offset  size  content
//   0       8     magic "CAMILFUL"
//   8       4     format version
//   12      4     key width in bits, 64
//   16      4     LCP level
//   20      4     nonzero when reverse complements were indexed
//   24      4     number of species
//   28      -     per species: 4 byte name length, then the name bytes
//   -       8     number of distinct cores
//   -       8     number of (core, species) pairs
//   -       4     ensemble bits
//   -       -     per sub table: 4 byte bits, 4 byte count, 4 byte present
//                 flag, then the used bitmap and the bucket array when present
//   -       -     the species array
// ---------------------------------------------------------------------------

#ifndef CAMIL_CFULL_H
#define CAMIL_CFULL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "logger.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define CFULL_MAGIC   "CAMILFUL"
#define CFULL_VERSION 1u

// split of the packed map value
#define CFULL_OFFSET_BITS 40
#define CFULL_COUNT_BITS  24
#define CFULL_MAX_PAIRS   (((uint64_t)1 << CFULL_OFFSET_BITS) - 1)
#define CFULL_MAX_COUNT   (((uint32_t)1 << CFULL_COUNT_BITS) - 1)

// entries one sub table is allowed to hold before the ensemble is split
// further. Keeps each sub table's resize cheap and its bucket count far
// below the 32 bit limit.
#define CFULL_SUB_TARGET ((uint64_t)1 << 26)
#define CFULL_MAX_ENSEMBLE_BITS 16

static inline uint64_t cfull_pack(uint64_t offset, uint32_t count) {
	return offset | ((uint64_t)count << CFULL_OFFSET_BITS);
}

static inline uint64_t cfull_offset(uint64_t value) {
	return value & CFULL_MAX_PAIRS;
}

static inline uint32_t cfull_count(uint64_t value) {
	return (uint32_t)(value >> CFULL_OFFSET_BITS);
}

// defined by the khashl instantiation in cfull.c
struct cfull_map;
struct camil_index;

struct cfull {
	int lcp_level;         // LCP level every core was parsed at
	int use_rc;            // nonzero when reverse complements were indexed
	uint32_t ngenomes;     // number of species
	char **names;          // display name per species, indexed by species id
	struct cfull_map *map; // key -> packed (offset, count)
	camil_sid *sids;       // the flat species array
	uint64_t nsids;        // pairs stored in sids
	uint64_t ncores;       // distinct cores the map is sized for
};

// Sizes the map for `ncores` distinct cores and the species array for
// `npairs` entries. The header fields are left as they are. Returns 0 on
// success, -1 when the counts exceed the format or memory runs out.
int cfull_alloc(struct cfull *full, uint64_t ncores, uint64_t npairs);

// Records that core `key` occupies sids[offset .. offset+count). Returns 0 on
// success, -1 on a repeated key or when out of memory.
int cfull_put(struct cfull *full, uint64_t key, uint64_t offset, uint32_t count);

// Looks a core up. Returns 0 and fills offset and count when present, -1 when
// absent.
int cfull_get(const struct cfull *full, uint64_t key, uint64_t *offset, uint32_t *count);

// Number of distinct cores stored.
uint64_t cfull_size(const struct cfull *full);

// Approximate heap footprint in bytes.
uint64_t cfull_memory(const struct cfull *full);

// Releases everything, names included.
void cfull_destroy(struct cfull *full);

// Writes the index in the format documented above. Returns 0 on success.
int cfull_save(const struct cfull *full, const char *path);

// Reads an index written by cfull_save(). Returns 0 on success.
int cfull_load(struct cfull *full, const char *path);

// Name of species sid, or "?" when the id is out of range.
const char *cfull_genome_name(const struct cfull *full, camil_sid sid);

// Writes a short human readable description of the index to stderr.
void cfull_report(const struct cfull *full);

// Builds a classify index over the named species. The output numbers them
// 0..nnames-1 in the order given. A core is kept when it occurs in at least
// one of them and in at most `max_share` of them; species outside the subset
// are ignored entirely. Fails on an unknown or repeated name. `threads`
// reducers sort the result. Returns 0 on success.
int cfull_subset(const struct cfull *full, char *const *names, uint32_t nnames, uint32_t max_share, int threads, struct camil_index *index);

#ifdef __cplusplus
}
#endif

#endif
