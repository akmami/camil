#include "camil.h"
#include "cfull.h"
#include "classify.h"
#include "index.h"
#include "logger.h"
#include "opt.h"
#include "seqio.h"
#include <stdio.h>
#include <stdlib.h>


static int check_inputs(char *const *paths, uint32_t count, const char *kind) {
	uint32_t i;
	int ok = 1;

	for (i = 0; i < count; i++) {
		if (!seq_readable(paths[i])) {
			log_error("cannot read the %s file %s", kind, paths[i]);
			ok = 0;
		}
	}
	return ok ? 0 : -1;
}

static FILE *open_output(const char *path, FILE *fallback) {
	FILE *out;

	if (path == NULL) {
		return fallback;
	}
	out = fopen(path, "w");
	if (out == NULL) {
		log_error("cannot write to %s", path);
	}
	return out;
}

static int classify_reads(const struct camil_index *index, const struct camil_opts *opts) {
	FILE *per_read = NULL;
	FILE *summary;
	struct camil_stats total;
	uint32_t i;
	int status = 0;

	if (opts->per_read_out != NULL) {
		per_read = open_output(opts->per_read_out, NULL);
		if (per_read == NULL) {
			return -1;
		}
		camil_classify_write_header(per_read);
	}

	summary = open_output(opts->summary_out, stdout);
	if (summary == NULL) {
		if (per_read != NULL) {
			fclose(per_read);
		}
		return -1;
	}

	if (camil_stats_init(&total, index->ngenomes) != 0) {
		status = -1;
		goto done;
	}

	for (i = 0; i < opts->nreads; i++) {
		struct camil_stats stats;

		if (camil_stats_init(&stats, index->ngenomes) != 0) {
			status = -1;
			break;
		}

		if (camil_classify_file(index, opts->reads[i], &opts->thresholds, opts->threads,
		                        per_read, &stats) != 0) {
			log_warn("classification of %s failed, continuing with the next file",
			         opts->reads[i]);
			camil_stats_destroy(&stats);
			status = -1;
			continue;
		}

		log_info("classified %s: %llu reads", opts->reads[i],
		         (unsigned long long)stats.total);
		camil_stats_report(&stats, index, opts->reads[i], summary);
		camil_stats_merge(&total, &stats);
		camil_stats_destroy(&stats);
	}

	// A combined table only carries information when several files were read.
	if (opts->nreads > 1) {
		camil_stats_report(&total, index, "all read files", summary);
	}
	camil_stats_destroy(&total);

done:
	if (per_read != NULL) {
		fclose(per_read);
	}
	if (summary != stdout) {
		fclose(summary);
	}
	return status;
}

// camil index: build the full index of the reference genomes and save it
static int command_index(const struct camil_opts *opts) {
	struct cfull full;

	// Genomes are not probed up front: with tens of thousands of them that is
	// a second pass over the file system, and a file that cannot be opened is
	// reported and skipped when its turn comes.
	if (camil_full_build(&full, opts->genomes, opts->ngenomes, opts->lcp_level, opts->use_rc, opts->threads) != 0) {
		return -1;
	}

	cfull_report(&full);

	if (cfull_save(&full, opts->index_out) != 0) {
		cfull_destroy(&full);
		return -1;
	}

	cfull_destroy(&full);
	return 0;
}

// camil subset: load a full index and write a classify index over some of its species
static int command_subset(const struct camil_opts *opts) {
	struct cfull full;
	struct camil_index index;

	if (cfull_load(&full, opts->index_in) != 0) {
		return -1;
	}
	cfull_report(&full);

	if (cfull_subset(&full, opts->species, opts->nspecies, (uint32_t)opts->max_share, opts->threads, &index) != 0) {
		cfull_destroy(&full);
		return -1;
	}
	cfull_destroy(&full);

	camil_index_report(&index);

	if (camil_index_save(&index, opts->index_out) != 0) {
		camil_index_destroy(&index);
		return -1;
	}

	camil_index_destroy(&index);
	return 0;
}

// camil classify: load an index and classify the given read files
static int command_classify(const struct camil_opts *opts) {
	struct camil_index index;
	int status;

	if (check_inputs(opts->reads, opts->nreads, "read") != 0) {
		return -1;
	}

	if (camil_index_load(&index, opts->index_in) != 0) {
		return -1;
	}

	camil_index_report(&index);
	status = classify_reads(&index, opts);
	camil_index_destroy(&index);
	return status;
}

// camil run: build an index in memory and classify in the same process
static int command_run(const struct camil_opts *opts) {
	struct camil_index index;
	int status;

	if (check_inputs(opts->reads, opts->nreads, "read") != 0) {
		return -1;
	}

	if (camil_index_build(&index, opts->genomes, opts->ngenomes, opts->lcp_level, opts->max_share, opts->use_rc, opts->threads) != 0) {
		return -1;
	}

	camil_index_report(&index);

	if (opts->index_out != NULL && camil_index_save(&index, opts->index_out) != 0) {
		camil_index_destroy(&index);
		return -1;
	}

	status = classify_reads(&index, opts);
	camil_index_destroy(&index);
	return status;
}

int main(int argc, char **argv) {
	struct camil_opts opts;
	int parsed;
	int status;

	parsed = camil_opts_parse(&opts, argc, argv);
	if (parsed != 0) {
		camil_opts_free(&opts);
		return parsed > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
	}

	log_set_verbose(opts.verbose);
	log_info(CAMIL_NAME " " CAMIL_VERSION " starting");

	LCP_INIT();

	switch (opts.command) {
	case CAMIL_CMD_INDEX:
		status = command_index(&opts);
		break;
	case CAMIL_CMD_CLASSIFY:
		status = command_classify(&opts);
		break;
	case CAMIL_CMD_RUN:
		status = command_run(&opts);
		break;
	case CAMIL_CMD_SUBSET:
		status = command_subset(&opts);
		break;
	default:
		camil_usage(stderr);
		status = -1;
		break;
	}

	camil_opts_free(&opts);

	log_info("done in %.3f seconds", log_elapsed());
	return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
