// ---------------------------------------------------------------------------
// index.h -- construction, storage and serialization of a camil core index.
//
// An index is the discriminative summary of a set of reference genomes:
//
//   1. Collection. Every reference is read sequence by sequence and parsed into
//      LCP cores at the requested level. Each core is reduced to a 64 bit key,
//      the lcptools label in the high half and a hash of the core's own bases
//      in the low half, and recorded against the species it came from. Unless
//      reverse complements are disabled each sequence is parsed twice, so that
//      reads from either strand match.
//
//   2. Reduction. Records are sorted, repeats of a core within one species
//      collapse to a single entry, and cores present in more than `max_share`
//      species are dropped. With the default of 1 only cores unique to a single
//      reference survive. Pruning only helps when references share cores, so an
//      unrelated reference set keeps nearly everything.
//
// Sequences are parsed in parallel and the reduction is parallel over disjoint
// key ranges, so neither phase has a lock on its hot path. The amount of
// sequence data held in flight is capped independently of the thread count so
// that a handful of large chromosomes cannot exhaust memory.
//
// On disk format
// --------------
// Values are written in host byte order; a magic number, a version and a key
// width let the loader reject anything an incompatible build produced.
//
//   offset  size  content
//   0       8     magic "CAMILIDX"
//   8       4     format version
//   12      4     key width in bits
//   16      4     LCP level the cores were computed at
//   20      4     max_share used when reducing
//   24      4     nonzero when reverse complements were indexed
//   28      4     margin folded into each key
//   32      4     number of species
//   36      -     per species: 4 byte name length, then the name bytes
//   -       4     directory bits
//   -       8     number of entries
//   -       -     directory, then keys, then values
//
// The three arrays are stored exactly as they sit in memory, so loading an
// index is a read and nothing else.
// ---------------------------------------------------------------------------

#ifndef CAMIL_INDEX_H
#define CAMIL_INDEX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "cbuild.h"
#include "ckey.h"
#include "ctable.h"
#include "logger.h"
#include "seqio.h"
#include "tpool.h"
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define CAMIL_INDEX_MAGIC   "CAMILIDX"
#define CAMIL_INDEX_VERSION 2u

// upper bound on the number of sequence bases queued for parsing at any one time
#define CAMIL_INDEX_INFLIGHT_BASES (1024ull * 1024ull * 1024ull)

// fraction of the key space an index may fill before accidental core matches
// are worth reporting
#define CAMIL_COLLISION_WARN 0.01

// representative core count per read, used only to turn the per core collision
// rate into a per read one in the warning
#define CAMIL_COLLISION_READ_CORES 20


struct camil_index {
	int lcp_level;       // LCP level used for every reference and every read
	int max_share;       // species a core may occur in and still be kept
	int use_rc;          // nonzero when reverse complements were indexed
	uint32_t margin;     // flanking bases folded into every key
	uint32_t ngenomes;   // number of references
	char **names;        // display name per species, indexed by species id
	struct ctable table; // frozen core key -> species, read only
};

// builds an index from ngenomes reference files. Each genome is reported under
// its short name when it has one, otherwise under its file name stripped of
// directory and extension.
int camil_index_build(struct camil_index *index, const struct camil_genome *genomes, uint32_t ngenomes, int lcp_level, int max_share, int use_rc, uint32_t margin, int threads);

// writes the index to path in the format documented above
int camil_index_save(const struct camil_index *index, const char *path);

// reads an index previously written by camil_index_save()
int camil_index_load(struct camil_index *index, const char *path);

// releases every resource held by the index
void camil_index_destroy(struct camil_index *index);

// name of species sid, or "?" when the id is out of range
const char *camil_index_genome_name(const struct camil_index *index, camil_sid sid);

// number of cores retained after reduction
uint64_t camil_index_size(const struct camil_index *index);

// writes a short human readable description of the index to stderr
void camil_index_report(const struct camil_index *index);

// warns when the index fills enough of the key space that chance matches start
// to matter. Called by camil_index_report()
void camil_index_warn_collisions(const struct camil_index *index);

#ifdef __cplusplus
}
#endif

#endif
