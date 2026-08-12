#ifndef CAMIL_CLASSIFY_H
#define CAMIL_CLASSIFY_H


#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "index.h"
#include "ctable.h"
#include "logger.h"
#include "seqio.h"
#include "tpool.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

// reads gathered before a batch is handed to the workers
#define CAMIL_BATCH_READS 8192u


// decision thresholds, see the discussion above
struct camil_thresholds {
	uint32_t min_hits; // minimum votes for the winning genome, default 1
	double min_ratio;  // minimum share of matched cores, default 1.0
};

// assigned has one counter per genome and is allocated by camil_stats_init()
struct camil_stats {
	uint64_t *assigned;    // reads attributed to each genome
	uint32_t ngenomes;     // length of assigned
	uint64_t ambiguous;    // matched something, but not decisively
	uint64_t unclassified; // no core of the read is in the index
	uint64_t total;        // reads seen
};

// fills thresholds with the defaults: one hit, full agreement
void camil_thresholds_default(struct camil_thresholds *thresholds);

// allocates the per genome counters
int camil_stats_init(struct camil_stats *stats, uint32_t ngenomes);

// releases the counters
void camil_stats_destroy(struct camil_stats *stats);

// adds every counter of src into dst. Used to fold per file tallies into a combined total
void camil_stats_merge(struct camil_stats *dst, const struct camil_stats *src);

// classifies every read in path against index
int camil_classify_file(const struct camil_index *index, const char *path, const struct camil_thresholds *thresholds, int threads, FILE *per_read, struct camil_stats *stats);

// writes the header line of the per read report
void camil_classify_write_header(FILE *per_read);

// writes a tally as a tab separated table: one row per genome, then the ambiguous, unclassified and total rows
void camil_stats_report(const struct camil_stats *stats, const struct camil_index *index, const char *label, FILE *out);

#ifdef __cplusplus
}
#endif

#endif
