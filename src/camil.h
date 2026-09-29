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

// species identifiers. Widening this to uint32_t is a one line change plus a
// bump of the index format version; nothing packs it into a fixed array.
typedef uint16_t camil_sid;

#define CAMIL_MAX_GENOMES 65534u
#define CAMIL_SID_NONE    ((camil_sid)0xFFFFu)

// a core may be shared by any number of species and still be kept; the limit is
// a run length now rather than the width of a field
#define CAMIL_MAX_SHARE CAMIL_MAX_GENOMES

// longest genome label a report can carry, terminator included
#define CAMIL_NAME_MAX 256

// one reference genome as it was requested on the command line
struct camil_genome {
	char *path;
	char *name;
};

// what the index stores against a core key. Four bytes, held in an array
// parallel to the keys so that a lookup that misses never touches it.
struct camil_value {
	camil_sid sid; // species the core was seen in
	uint8_t level; // LCP level the core was parsed at, for later multi-level use
	uint8_t flags; // reserved
};

// fails the build if the value record grows, which would silently undo the
// layout the lookup path relies on
typedef char camil_value_size_check[sizeof(struct camil_value) == 4 ? 1 : -1];

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

// finalizer used to mix a 64 bit accumulator, the splitmix64 tail
static inline uint64_t camil_mix64(uint64_t x) {
	x ^= x >> 33;
	x *= 0xFF51AFD7ED558CCDULL;
	x ^= x >> 33;
	x *= 0xC4CEB9FE1A85EC53ULL;
	x ^= x >> 33;
	return x;
}

#ifdef __cplusplus
}
#endif

#endif
