#ifndef CAMIL_OPT_H
#define CAMIL_OPT_H


#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "classify.h"
#include "logger.h"
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


// default number of worker threads
#define CAMIL_DEFAULT_THREADS 1

// default LCP level
#define CAMIL_DEFAULT_LEVEL 4

enum camil_command {
	CAMIL_CMD_NONE = 0,
	CAMIL_CMD_INDEX,
	CAMIL_CMD_CLASSIFY,
	CAMIL_CMD_RUN,
	CAMIL_CMD_HELP,
	CAMIL_CMD_VERSION
};

// everything the command line can express
struct camil_opts {
	enum camil_command command;

	int lcp_level; // -l, LCP level for references and reads alike
	int max_share; // -n, genomes a core may occur in and still be indexed
	int threads;   // -t, worker threads
	uint32_t margin; // --margin, flanking bases folded into every core key
	int use_rc;    // cleared by --no-rc
	int verbose;   // -v

	struct camil_thresholds thresholds; // --min-hits, --min-ratio

	const char *index_out;        // index: -o, run: --save-index
	const char *index_in;         // classify: -i
	const char *per_read_out;     // classify and run: -o
	const char *summary_out;      // classify and run: -s, stdout when NULL

	struct camil_genome *genomes; // index: positional, run: -g, either: -G
	uint32_t ngenomes;
	char **reads;                 // classify and run: positional
	uint32_t nreads;
};

// parses `argc`/`argv` into `opts`, filling in the defaults first.
int camil_opts_parse(struct camil_opts *opts, int argc, char **argv);

// releases the vectors owned by `opts`.
void camil_opts_free(struct camil_opts *opts);

// prints the top level usage summary.
void camil_usage(FILE *out);

// prints the usage of one subcommand, including its options and defaults.
void camil_usage_command(FILE *out, enum camil_command command);

#ifdef __cplusplus
}
#endif

#endif
