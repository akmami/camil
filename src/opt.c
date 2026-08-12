#include "opt.h"


// appends value to a growable vector of argv pointers, returns 0 on success.
static int vector_push(char ***vector, uint32_t *count, char *value) {
	char **grown = (char **)realloc(*vector, ((size_t)(*count) + 1) * sizeof(char *));

	if (grown == NULL) {
		log_error("out of memory while collecting arguments");
		return -1;
	}
	*vector = grown;
	(*vector)[*count] = value;
	(*count)++;
	return 0;
}

// matches arg against a short and a long spelling
static int opt_is(const char *arg, const char *short_form, const char *long_form) {
	return (short_form != NULL && strcmp(arg, short_form) == 0) || (long_form != NULL && strcmp(arg, long_form) == 0);
}

// consumes the value that follows an option
static char *opt_value(int argc, char **argv, int *i, const char *name) {
	if (*i + 1 >= argc) {
		log_error("option %s expects a value", name);
		return NULL;
	}
	(*i)++;
	return argv[*i];
}

// parses an integer
static int opt_int(const char *text, const char *name, long *out) {
	char *end = NULL;
	long value;

	errno = 0;
	value = strtol(text, &end, 10);

	if (end == text || *end != '\0') {
		log_error("option %s expects an integer, got '%s'", name, text);
		return -1;
	}
	if (errno == ERANGE || value < INT_MIN || value > INT_MAX) {
		log_error("option %s is out of range: '%s'", name, text);
		return -1;
	}
	*out = value;
	return 0;
}

// parses a floating point value in the same spirit as opt_int()
static int opt_double(const char *text, const char *name, double *out) {
	char *end = NULL;
	double value;

	errno = 0;
	value = strtod(text, &end);

	if (end == text || *end != '\0' || !(value == value)) {
		log_error("option %s expects a number, got '%s'", name, text);
		return -1;
	}
	if (errno == ERANGE) {
		log_error("option %s is out of range: '%s'", name, text);
		return -1;
	}
	*out = value;
	return 0;
}

static enum camil_command opt_command(const char *word) {
	if (strcmp(word, "index") == 0) {
		return CAMIL_CMD_INDEX;
	}
	if (strcmp(word, "classify") == 0) {
		return CAMIL_CMD_CLASSIFY;
	}
	if (strcmp(word, "run") == 0) {
		return CAMIL_CMD_RUN;
	}
	if (strcmp(word, "help") == 0 || strcmp(word, "-h") == 0 || strcmp(word, "--help") == 0) {
		return CAMIL_CMD_HELP;
	}
	if (strcmp(word, "version") == 0 || strcmp(word, "-V") == 0 ||
	    strcmp(word, "--version") == 0) {
		return CAMIL_CMD_VERSION;
	}
	return CAMIL_CMD_NONE;
}

void camil_usage(FILE *out) {
	fprintf(out,
	        CAMIL_NAME " (" CAMIL_VERSION
	        "): read classification with locally consistent parsing cores\n"
	        "\n"
	        "usage: ./" CAMIL_NAME " <command> [options]\n"
	        "\n"
	        "commands:\n"
	        "  index      build a core index from reference genomes\n"
	        "  classify   assign reads using a previously built index\n"
	        "  run        build an index in memory and classify reads in one go\n"
	        "  help       show this message, or 'help <command>' for details\n"
	        "  version    print the version and exit\n");
}

void camil_usage_command(FILE *out, enum camil_command command) {
	switch (command) {
	case CAMIL_CMD_INDEX:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " index -o <index> [options] <genome.fa> [genome.fa ...]\n"
		        "\n"
		        "Parses every reference genome into LCP cores and stores the cores that\n"
		        "occur in at most --max-share genomes. Genome names are derived from the\n"
		        "file names.\n"
		        "\n"
		        "options:\n"
		        "  -o, --output FILE     where to write the index (required)\n"
		        "  -l, --level INT       LCP level (default %d)\n"
		        "  -n, --max-share INT   genomes a core may occur in, 1 to %d (default 1)\n"
		        "  -t, --threads INT     worker threads (default %d)\n"
		        "      --no-rc           do not index reverse complements\n"
		        "  -v, --verbose         print debug messages\n"
		        "  -h, --help            show this message\n",
		        CAMIL_DEFAULT_LEVEL, CAMIL_MAX_SHARE, CAMIL_DEFAULT_THREADS);
		break;
	case CAMIL_CMD_CLASSIFY:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " classify -i <index> [options] <reads.fq> [reads.fq ...]\n"
		        "\n"
		        "Assigns each read to the genome that owns most of its cores. Reads with no\n"
		        "core in the index are unclassified; reads whose cores disagree, or that\n"
		        "fail a threshold, are ambiguous.\n"
		        "\n"
		        "options:\n"
		        "  -i, --index FILE      index written by '" CAMIL_NAME " index' (required)\n"
		        "  -o, --output FILE     per read report, tab separated\n"
		        "  -s, --summary FILE    summary table (default stdout)\n"
		        "  -t, --threads INT     worker threads (default %d)\n"
		        "      --min-hits INT    minimum cores supporting the winner (default 1)\n"
		        "      --min-ratio FLOAT share of matched cores the winner must hold,\n"
		        "                        between 0 and 1 (default 1.0)\n"
		        "  -v, --verbose         print debug messages\n"
		        "  -h, --help            show this message\n",
		        CAMIL_DEFAULT_THREADS);
		break;
	case CAMIL_CMD_RUN:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " run -g <genome.fa> [-g ...] [options] <reads.fq> [reads.fq ...]\n"
		        "\n"
		        "Builds an index in memory and immediately classifies the given reads.\n"
		        "Equivalent to 'index' followed by 'classify' without writing the index,\n"
		        "unless --save-index is given.\n"
		        "\n"
		        "options:\n"
		        "  -g, --genome FILE     reference genome, repeat once per genome (required)\n"
		        "  -o, --output FILE     per read report, tab separated\n"
		        "  -s, --summary FILE    summary table (default stdout)\n"
		        "      --save-index FILE also write the index to this file\n"
		        "  -l, --level INT       LCP level (default %d)\n"
		        "  -n, --max-share INT   genomes a core may occur in, 1 to %d (default 1)\n"
		        "  -t, --threads INT     worker threads (default %d)\n"
		        "      --no-rc           do not index reverse complements\n"
		        "      --min-hits INT    minimum cores supporting the winner (default 1)\n"
		        "      --min-ratio FLOAT share of matched cores the winner must hold,\n"
		        "                        between 0 and 1 (default 1.0)\n"
		        "  -v, --verbose         print debug messages\n"
		        "  -h, --help            show this message\n",
		        CAMIL_DEFAULT_LEVEL, CAMIL_MAX_SHARE, CAMIL_DEFAULT_THREADS);
		break;
	default:
		camil_usage(out);
		break;
	}
}

// fills opts with the defaults shared by every subcommand
static void opt_defaults(struct camil_opts *opts) {
	memset(opts, 0, sizeof(*opts));
	opts->lcp_level = CAMIL_DEFAULT_LEVEL;
	opts->max_share = 1;
	opts->threads = CAMIL_DEFAULT_THREADS;
	opts->use_rc = 1;
	camil_thresholds_default(&opts->thresholds);
}

// validates the option values and the presence of the arguments each subcommand requires
static int opt_validate(const struct camil_opts *opts) {
	if (opts->lcp_level < 1) {
		log_error("the LCP level must be at least 1");
		return -1;
	}
	if (opts->max_share < 1 || opts->max_share > CAMIL_MAX_SHARE) {
		log_error("--max-share must be between 1 and %d", CAMIL_MAX_SHARE);
		return -1;
	}
	if (opts->threads < 1) {
		log_error("--threads must be at least 1");
		return -1;
	}
	if (!(opts->thresholds.min_ratio >= 0.0 && opts->thresholds.min_ratio <= 1.0)) {
		log_error("--min-ratio must be between 0 and 1");
		return -1;
	}
	if (opts->thresholds.min_hits < 1) {
		log_error("--min-hits must be at least 1");
		return -1;
	}

	switch (opts->command) {
	case CAMIL_CMD_INDEX:
		if (opts->index_out == NULL) {
			log_error("index: --output is required");
			return -1;
		}
		if (opts->ngenomes == 0) {
			log_error("index: at least one reference genome is required");
			return -1;
		}
		break;
	case CAMIL_CMD_CLASSIFY:
		if (opts->index_in == NULL) {
			log_error("classify: --index is required");
			return -1;
		}
		if (opts->nreads == 0) {
			log_error("classify: at least one read file is required");
			return -1;
		}
		break;
	case CAMIL_CMD_RUN:
		if (opts->ngenomes == 0) {
			log_error("run: at least one --genome is required");
			return -1;
		}
		if (opts->nreads == 0) {
			log_error("run: at least one read file is required");
			return -1;
		}
		break;
	default:
		break;
	}

	return 0;
}

int camil_opts_parse(struct camil_opts *opts, int argc, char **argv) {
	int i;
	int only_positional = 0;

	opt_defaults(opts);

	if (argc < 2) {
		camil_usage(stderr);
		return -1;
	}

	opts->command = opt_command(argv[1]);
	if (opts->command == CAMIL_CMD_NONE) {
		log_error("unknown command '%s'", argv[1]);
		camil_usage(stderr);
		return -1;
	}

	if (opts->command == CAMIL_CMD_VERSION) {
		printf(CAMIL_NAME " " CAMIL_VERSION "\n");
		return 1;
	}

	if (opts->command == CAMIL_CMD_HELP) {
		// 'help <command>' prints the details of that command.
		if (argc > 2) {
			camil_usage_command(stdout, opt_command(argv[2]));
		} else {
			camil_usage(stdout);
		}
		return 1;
	}

	for (i = 2; i < argc; i++) {
		char *arg = argv[i];
		char *value;
		long number;
		double ratio;

		if (only_positional || arg[0] != '-' || arg[1] == '\0') {
			// positional: a genome for 'index', a read file otherwise.
			if (opts->command == CAMIL_CMD_INDEX) {
				if (vector_push(&opts->genomes, &opts->ngenomes, arg) != 0) {
					return -1;
				}
			} else if (vector_push(&opts->reads, &opts->nreads, arg) != 0) {
				return -1;
			}
			continue;
		}

		if (strcmp(arg, "--") == 0) {
			only_positional = 1;
			continue;
		}

		if (opt_is(arg, "-h", "--help")) {
			camil_usage_command(stdout, opts->command);
			return 1;
		}
		if (opt_is(arg, "-v", "--verbose")) {
			opts->verbose = 1;
			continue;
		}
		if (opt_is(arg, NULL, "--no-rc")) {
			opts->use_rc = 0;
			continue;
		}

		if (opt_is(arg, "-l", "--level")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    opt_int(value, arg, &number) != 0) {
				return -1;
			}
			opts->lcp_level = (int)number;
			continue;
		}
		if (opt_is(arg, "-n", "--max-share")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    opt_int(value, arg, &number) != 0) {
				return -1;
			}
			opts->max_share = (int)number;
			continue;
		}
		if (opt_is(arg, "-t", "--threads")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    opt_int(value, arg, &number) != 0) {
				return -1;
			}
			opts->threads = (int)number;
			continue;
		}
		if (opt_is(arg, NULL, "--min-hits")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    opt_int(value, arg, &number) != 0) {
				return -1;
			}
			opts->thresholds.min_hits = number < 0 ? 0 : (uint32_t)number;
			continue;
		}
		if (opt_is(arg, NULL, "--min-ratio")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    opt_double(value, arg, &ratio) != 0) {
				return -1;
			}
			opts->thresholds.min_ratio = ratio;
			continue;
		}
		if (opt_is(arg, "-i", "--index")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			opts->index_in = value;
			continue;
		}
		if (opt_is(arg, "-s", "--summary")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			opts->summary_out = value;
			continue;
		}
		if (opt_is(arg, NULL, "--save-index")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			opts->index_out = value;
			continue;
		}
		if (opt_is(arg, "-g", "--genome")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			if (vector_push(&opts->genomes, &opts->ngenomes, value) != 0) {
				return -1;
			}
			continue;
		}
		if (opt_is(arg, "-o", "--output")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			// For 'index' the output is the index itself; for the other
			// commands it is the per read report.
			if (opts->command == CAMIL_CMD_INDEX) {
				opts->index_out = value;
			} else {
				opts->per_read_out = value;
			}
			continue;
		}

		log_error("unknown option '%s'", arg);
		camil_usage_command(stderr, opts->command);
		return -1;
	}

	if (opt_validate(opts) != 0) {
		return -1;
	}

	return 0;
}

void camil_opts_free(struct camil_opts *opts) {
	free(opts->genomes);
	free(opts->reads);
	opts->genomes = NULL;
	opts->reads = NULL;
	opts->ngenomes = 0;
	opts->nreads = 0;
}
