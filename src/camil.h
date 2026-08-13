#ifndef CAMIL_H
#define CAMIL_H

#ifdef __cplusplus
extern "C" {
#endif


#include "core.h"
#include "lps.h"
#include "encoding.h"
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


#define CAMIL_NAME    "camil"
#define CAMIL_VERSION "1.0.0"

// largest number of genomes a single core may belong to and still be kept
#define CAMIL_MAX_SHARE 3

// genome identifiers are single bytes, which caps a run at 255 references
typedef uint8_t camil_gid;

#define CAMIL_MAX_GENOMES 255u

// sentinel stored in ngenomes while indexing, marking a core that turned up in more than CAMIL_MAX_SHARE genomes
#define CAMIL_OVERFLOW 255u

// longest genome label a report can carry, terminator included
#define CAMIL_NAME_MAX 256

// one reference genome as it was requested on the command line
struct camil_genome {
	char *path;
	char *name;
};

// one core of the index
struct camil_entry {
	lcp_label label;                 // LCP core label, the key
	uint8_t ngenomes;                // 0 empty, 1..CAMIL_MAX_SHARE, or CAMIL_OVERFLOW
	camil_gid gids[CAMIL_MAX_SHARE]; // genomes the core occurs in
};

// fails the build if the entry ever grows past a cache line quarter, which would silently undo the layout the lookup path relies on
typedef char camil_entry_size_check[sizeof(struct camil_entry) <= 16 ? 1 : -1];

// cache line size assumed when padding per thread scratch space. Being wrong here costs a little memory, never correctness.
#define CAMIL_CACHE_LINE 64

// outcome of classifying a single read
enum camil_status {
	CAMIL_UNCLASSIFIED = 0, 	// none of the read cores are present in the index
	CAMIL_AMBIGUOUS    = 1,		// the read has hits but they fail the thresholds, or two or more genomes are tied for the best score
	CAMIL_ASSIGNED     = 2		// the read was attributed to exactly one genome
};

// human readable form of enum camil_status, never NULL
const char *camil_status_name(int status);

// duplicates a NUL terminated string, NULL when out of memory.
//
// strdup is POSIX rather than ISO C, so a strict -std=cNN hides it and the
// call silently becomes an implicit int, truncating the pointer. Carrying our
// own keeps the build independent of feature test macros and include order.
static inline char *camil_strdup(const char *text) {
	size_t size = strlen(text) + 1;
	char *copy = (char *)malloc(size);

	if (copy != NULL) {
		memcpy(copy, text, size);
	}
	return copy;
}

// mixes a core label into a 64 bit hash
static inline uint64_t camil_hash(lcp_label label) {
	uint64_t x = (uint64_t)label;

	x += 0x9E3779B97F4A7C15ULL;
	x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
	x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
	x = x ^ (x >> 31);
	return x;
}

#ifdef __cplusplus
}
#endif

#endif
