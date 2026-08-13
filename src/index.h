// ---------------------------------------------------------------------------
// index.h -- construction, storage and serialization of a camil core index.
//
// An index is the discriminative summary of a set of reference genomes. It is
// produced in two steps:
//
//   1. Collection. Every reference is read sequence by sequence and parsed
//      into LCP cores at the requested level. Each distinct core label is
//      recorded in the core table together with the genome it came from. A
//      label seen a thousand times in one genome is stored once; a label seen
//      in several genomes accumulates several genome ids.
//
//      Unless reverse complement handling is disabled, each reference sequence
//      is parsed twice, once as given and once as its reverse complement, so
//      that reads sequenced from either strand can match.
//
//   2. Pruning. Labels present in more than max_share genomes are dropped.
//      With the default max_share of 1 only labels unique to a single genome
//      survive, which reproduces the strict behaviour of the original script.
//      Raising it to 2 or 3 keeps cores shared by a few references, trading
//      specificity for sensitivity. CAMIL_MAX_SHARE is the hard ceiling.
//
// Sequences are parsed in parallel. Because parsing dominates the runtime and
// the core table is sharded, indexing scales close to linearly with the thread
// count until input decompression becomes the bottleneck. The amount of
// sequence data held in flight is capped independently of the thread count so
// that a handful of large chromosomes cannot exhaust memory.
//
// On disk format
// --------------
// Index files are little endian in practice: values are written with the host
// byte order, and a magic number plus a version and a type width field let the
// loader reject files produced by an incompatible build (for instance one
// compiled with 32 bit core labels reading a 64 bit index).
//
//   offset  size  content
//   0       8     magic "CAMILIDX"
//   8       4     format version
//   12      4     width of lcp_label in bits (32 or 64)
//   16      4     LCP level the cores were computed at
//   20      4     max_share used when pruning
//   24      4     nonzero when reverse complements were indexed
//   28      4     number of genomes
//   32      -     per genome: 4 byte name length, then the name bytes
//   -       8     table capacity in slots
//   -       8     number of occupied slots
//   -       -     the slot array, exactly as it sits in memory
//
// The table is stored raw rather than as a list of entries, which makes
// loading a single read with no rehashing: an index of a few hundred million
// cores comes back in the time it takes to stream it off disk.
// ---------------------------------------------------------------------------

#ifndef CAMIL_INDEX_H
#define CAMIL_INDEX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "ctable.h"
#include "cbuild.h"
#include "logger.h"
#include "seqio.h"
#include "tpool.h"
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define CAMIL_INDEX_MAGIC   "CAMILIDX"
#define CAMIL_INDEX_VERSION 1u

// upper bound on the number of sequence bases queued for parsing at any one time
#define CAMIL_INDEX_INFLIGHT_BASES (1024ull * 1024ull * 1024ull)


struct camil_index {
	int lcp_level;       // LCP level used for every reference and every read
	int max_share;       // genomes a core may occur in and still be kept
	int use_rc;          // nonzero when reverse complements were indexed
	uint32_t ngenomes;   // number of references
	char **names;        // display name per genome, indexed by genome id
	struct ctable table; // frozen core label -> genome ids, read only
};

// builds an index from ngenomes reference files
int camil_index_build(struct camil_index *index, const struct camil_genome *genomes, uint32_t ngenomes, int lcp_level, int max_share, int use_rc, int threads);

// writes the index to path in the format documented above
int camil_index_save(const struct camil_index *index, const char *path);

// reads an index previously written by camil_index_save()
int camil_index_load(struct camil_index *index, const char *path);

// releases every resource held by the index
void camil_index_destroy(struct camil_index *index);

// name of genome gid, or "?" when the id is out of range
const char *camil_index_genome_name(const struct camil_index *index, camil_gid gid);

// number of cores retained after pruning
uint64_t camil_index_size(const struct camil_index *index);

// writes a short human readable description of the index to stderr: level, sharing limit, genome names and the number of retained cores
void camil_index_report(const struct camil_index *index);

#ifdef __cplusplus
}
#endif

#endif
