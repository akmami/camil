#include "opt.h"


// longest line a --genome-list file may contain, newline included
#define CAMIL_LIST_LINE 8192

// duplicates a byte range as a NUL terminated string
static char *dup_range(const char *begin, size_t len) {
	char *copy = (char *)malloc(len + 1);

	if (copy == NULL) {
		log_error("out of memory while collecting arguments");
		return NULL;
	}
	memcpy(copy, begin, len);
	copy[len] = '\0';
	return copy;
}

// drops trailing blanks from a mutable string
static void trim_end(char *text) {
	size_t len = strlen(text);

	while (len > 0 && (text[len - 1] == ' ' || text[len - 1] == '\t')) {
		text[--len] = '\0';
	}
}

// appends a genome, taking ownership of both strings
static int genome_push(struct camil_opts *opts, char *path, char *name) {
	struct camil_genome *grown;

	grown = (struct camil_genome *)realloc(opts->genomes, ((size_t)opts->ngenomes + 1) * sizeof(struct camil_genome));
	if (grown == NULL) {
		log_error("out of memory while collecting genomes");
		free(path);
		free(name);
		return -1;
	}
	opts->genomes = grown;
	opts->genomes[opts->ngenomes].path = path;
	opts->genomes[opts->ngenomes].name = name;
	opts->ngenomes++;
	return 0;
}

// accepts "path" or "path,shortname"
static int genome_add_spec(struct camil_opts *opts, const char *spec) {
	const char *comma = strrchr(spec, ',');
	char *path;
	char *name = NULL;

	if (comma != NULL && comma != spec && comma[1] != '\0') {
		path = dup_range(spec, (size_t)(comma - spec));
		name = dup_range(comma + 1, strlen(comma + 1));
		if (path == NULL || name == NULL) {
			free(path);
			free(name);
			return -1;
		}
	} else {
		path = dup_range(spec, strlen(spec));
		if (path == NULL) {
			return -1;
		}
	}

	return genome_push(opts, path, name);
}

// reads a tab separated list of genomes: one path per line, optionally followed by a tab and the short name to report it under
static int genome_add_list(struct camil_opts *opts, const char *path) {
	char line[CAMIL_LIST_LINE];
	FILE *in;
	unsigned long lineno = 0;
	uint32_t before = opts->ngenomes;

	in = fopen(path, "r");
	if (in == NULL) {
		log_error("cannot read the genome list %s", path);
		return -1;
	}

	while (fgets(line, sizeof(line), in) != NULL) {
		char *begin = line;
		char *tab;
		char *name = NULL;
		char *path_copy;
		char *name_copy = NULL;
		size_t len = strlen(line);

		lineno++;

		// full buffer with no newline means the line was cut in half, and silently indexing half a path would be worse than stopping
		if (len + 1 == sizeof(line) && line[len - 1] != '\n') {
			log_error("%s: line %lu is longer than %d characters", path, lineno, (int)sizeof(line) - 1);
			fclose(in);
			return -1;
		}

		while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
			line[--len] = '\0';
		}
		while (*begin == ' ' || *begin == '\t') {
			begin++;
		}
		if (*begin == '\0' || *begin == '#') {
			continue;
		}

		tab = strchr(begin, '\t');
		if (tab != NULL) {
			*tab = '\0';
			name = tab + 1;
			while (*name == ' ' || *name == '\t') {
				name++;
			}
			trim_end(name);
			if (*name == '\0') {
				name = NULL;
			}
		}
		trim_end(begin);
		if (*begin == '\0') {
			continue;
		}

		path_copy = dup_range(begin, strlen(begin));
		if (name != NULL) {
			name_copy = dup_range(name, strlen(name));
		}
		if (path_copy == NULL || (name != NULL && name_copy == NULL)) {
			free(path_copy);
			free(name_copy);
			fclose(in);
			return -1;
		}
		if (genome_push(opts, path_copy, name_copy) != 0) {
			fclose(in);
			return -1;
		}
	}

	if (ferror(in)) {
		log_error("failed while reading the genome list %s", path);
		fclose(in);
		return -1;
	}
	fclose(in);

	log_info("read %u genome(s) from %s", opts->ngenomes - before, path);
	return 0;
}

// appends a species name, taking a copy
static int species_add(struct camil_opts *opts, const char *name) {
	char **grown;
	char *copy;

	copy = dup_range(name, strlen(name));
	if (copy == NULL) {
		return -1;
	}
	grown = (char **)realloc(opts->species, ((size_t)opts->nspecies + 1) * sizeof(char *));
	if (grown == NULL) {
		log_error("out of memory while collecting species names");
		free(copy);
		return -1;
	}
	opts->species = grown;
	opts->species[opts->nspecies++] = copy;
	return 0;
}

// reads species names from a file, one per line. Blank lines and lines
// starting with '#' are ignored; only the first tab separated field counts, so
// a genome list written for --genome-list can be reused when its names are
// in the first column.
static int species_add_list(struct camil_opts *opts, const char *path) {
	char line[CAMIL_LIST_LINE];
	FILE *in;
	unsigned long lineno = 0;
	uint32_t before = opts->nspecies;

	in = fopen(path, "r");
	if (in == NULL) {
		log_error("cannot read the species list %s", path);
		return -1;
	}

	while (fgets(line, sizeof(line), in) != NULL) {
		char *begin = line;
		char *tab;
		size_t len = strlen(line);

		lineno++;

		if (len + 1 == sizeof(line) && line[len - 1] != '\n') {
			log_error("%s: line %lu is longer than %d characters", path, lineno, (int)sizeof(line) - 1);
			fclose(in);
			return -1;
		}

		while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
			line[--len] = '\0';
		}
		while (*begin == ' ' || *begin == '\t') {
			begin++;
		}
		if (*begin == '\0' || *begin == '#') {
			continue;
		}
		tab = strchr(begin, '\t');
		if (tab != NULL) {
			*tab = '\0';
		}
		trim_end(begin);
		if (*begin == '\0') {
			continue;
		}
		if (species_add(opts, begin) != 0) {
			fclose(in);
			return -1;
		}
	}

	if (ferror(in)) {
		log_error("failed while reading the species list %s", path);
		fclose(in);
		return -1;
	}
	fclose(in);

	log_info("read %u species name(s) from %s", opts->nspecies - before, path);
	return 0;
}

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
	if (strcmp(word, "subset") == 0) {
		return CAMIL_CMD_SUBSET;
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
	        "usage: ./" CAMIL_NAME " <command> [options]\n"
	        "\n"
	        "commands:\n"
	        "  index      build the full core index of a set of reference genomes\n"
	        "  subset     extract a classify index for some of its species\n"
	        "  classify   assign reads using a classify index\n"
	        "  run        build a classify index in memory and classify reads in one go\n"
	        "  help       show this message, or 'help <command>' for details\n"
	        "  version    print the version and exit\n");
}

void camil_usage_command(FILE *out, enum camil_command command) {
	switch (command) {
	case CAMIL_CMD_INDEX:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " index -o <index> [options] <genome.fa[,name]> ...\n"
		        "\n"
		        "Parses every reference genome into LCP cores and stores every core with\n"
		        "the complete list of genomes it occurs in. Nothing is filtered: the\n"
		        "result is the full index that '" CAMIL_NAME " subset' carves classify\n"
		        "indexes out of, so the references are parsed once for any number of\n"
		        "subsets.\n"
		        "\n"
		        "A genome is reported under the short name written after a comma, or under\n"
		        "its file name with the directory and extension removed when no name is\n"
		        "given. '/path/to/GRCh38.fa.gz,human' is reported as 'human'.\n"
		        "\n"
		        "options:\n"
		        "  -o, --output FILE         where to write the full index (required)\n"
		        "  -G, --genome-list FILE    read genomes from a tab separated file, one\n"
		        "                            'path <tab> name' per line, the name optional\n"
		        "  -l, --level INT           LCP level (default %d)\n"
		        "  -t, --threads INT         worker threads (default %d)\n"
		        "      --no-rc               do not index reverse complements\n"
		        "  -v, --verbose             print debug messages\n"
		        "  -h, --help                show this message\n",
		        CAMIL_DEFAULT_LEVEL, CAMIL_DEFAULT_THREADS);
		break;
	case CAMIL_CMD_SUBSET:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " subset -i <full index> -o <index> [options] <name> ...\n"
		        "\n"
		        "Builds a classify index over some of the species of a full index. A core\n"
		        "is kept when it occurs in at least one of the named species and in at\n"
		        "most --max-share of them; species that were not named are ignored\n"
		        "entirely, so a core that is unique among the chosen species is kept even\n"
		        "when the full index shares it widely. The output numbers the species\n"
		        "0, 1, ... in the order they were named.\n"
		        "\n"
		        "options:\n"
		        "  -i, --index FILE          full index written by '" CAMIL_NAME " index' (required)\n"
		        "  -o, --output FILE         where to write the classify index (required)\n"
		        "  -S, --species-list FILE   read species names from a file, one per line;\n"
		        "                            only the first tab separated field is used\n"
		        "  -n, --max-share INT       chosen species a core may occur in (default 1)\n"
		        "  -t, --threads INT         worker threads (default %d)\n"
		        "  -v, --verbose             print debug messages\n"
		        "  -h, --help                show this message\n",
		        CAMIL_DEFAULT_THREADS);
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
		        "  -i, --index FILE          classify index written by '" CAMIL_NAME " subset' (required)\n"
		        "  -o, --output FILE         per read report, tab separated\n"
		        "  -s, --summary FILE        summary table (default stdout)\n"
		        "  -t, --threads INT         worker threads (default %d)\n"
		        "      --min-hits INT        minimum cores supporting the winner (default 1)\n"
		        "      --min-ratio FLOAT     share of matched cores the winner must hold,\n"
		        "                            between 0 and 1 (default 1.0)\n"
		        "  -v, --verbose             print debug messages\n"
		        "  -h, --help                show this message\n",
		        CAMIL_DEFAULT_THREADS);
		break;
	case CAMIL_CMD_RUN:
		fprintf(out,
		        "usage: ./" CAMIL_NAME " run -g <genome.fa[,name]> [-g ...] [options] <reads.fq> ...\n"
		        "\n"
		        "Builds an index in memory and immediately classifies the given reads.\n"
		        "Equivalent to 'index' followed by 'classify' without writing the index,\n"
		        "unless --save-index is given.\n"
		        "\n"
		        "A genome is reported under the short name written after a comma, or under\n"
		        "its file name with the directory and extension removed when no name is\n"
		        "given. '-g /data/refs/GRCh38.fa.gz,Human' is reported as 'Human'.\n"
		        "\n"
		        "options:\n"
		        "  -g, --genome FILE[,NAME]  reference genome, repeat once per genome\n"
		        "  -G, --genome-list FILE    read genomes from a tab separated file, one\n"
		        "                            'path <tab> name' per line, the name optional\n"
		        "  -o, --output FILE         per read report, tab separated\n"
		        "  -s, --summary FILE        summary table (default stdout)\n"
		        "      --save-index FILE     also write the index to this file\n"
		        "  -l, --level INT           LCP level (default %d)\n"
		        "  -n, --max-share INT       species a core may occur in (default 1)\n"
		        "  -t, --threads INT         worker threads (default %d)\n"
		        "      --no-rc               do not index reverse complements\n"
		        "      --min-hits INT        minimum cores supporting the winner (default 1)\n"
		        "      --min-ratio FLOAT     share of matched cores the winner must hold,\n"
		        "                            between 0 and 1 (default 1.0)\n"
		        "  -v, --verbose             print debug messages\n"
		        "  -h, --help                show this message\n",
		        CAMIL_DEFAULT_LEVEL, CAMIL_DEFAULT_THREADS);
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
	if (opts->max_share < 1 || (uint32_t)opts->max_share > CAMIL_MAX_SHARE) {
		log_error("--max-share must be between 1 and %u", CAMIL_MAX_SHARE);
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
	case CAMIL_CMD_SUBSET:
		if (opts->index_in == NULL) {
			log_error("subset: --index is required");
			return -1;
		}
		if (opts->index_out == NULL) {
			log_error("subset: --output is required");
			return -1;
		}
		if (opts->nspecies == 0) {
			log_error("subset: at least one species name is required");
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
			// positional: a genome for 'index', a species name for 'subset', a read file otherwise.
			if (opts->command == CAMIL_CMD_INDEX) {
				if (genome_add_spec(opts, arg) != 0) {
					return -1;
				}
			} else if (opts->command == CAMIL_CMD_SUBSET) {
				if (species_add(opts, arg) != 0) {
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
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    genome_add_spec(opts, value) != 0) {
				return -1;
			}
			continue;
		}
		if (opt_is(arg, "-S", "--species-list")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    species_add_list(opts, value) != 0) {
				return -1;
			}
			continue;
		}
		if (opt_is(arg, "-G", "--genome-list")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL ||
			    genome_add_list(opts, value) != 0) {
				return -1;
			}
			continue;
		}
		if (opt_is(arg, "-o", "--output")) {
			if ((value = opt_value(argc, argv, &i, arg)) == NULL) {
				return -1;
			}
			// For 'index' and 'subset' the output is the index itself; for
			// the other commands it is the per read report.
			if (opts->command == CAMIL_CMD_INDEX || opts->command == CAMIL_CMD_SUBSET) {
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
	uint32_t i;

	// genome paths and names are copies, unlike the read files, which point straight into argv
	for (i = 0; i < opts->ngenomes; i++) {
		free(opts->genomes[i].path);
		free(opts->genomes[i].name);
	}
	for (i = 0; i < opts->nspecies; i++) {
		free(opts->species[i]);
	}
	free(opts->genomes);
	free(opts->reads);
	free(opts->species);
	opts->genomes = NULL;
	opts->reads = NULL;
	opts->species = NULL;
	opts->ngenomes = 0;
	opts->nreads = 0;
	opts->nspecies = 0;
}
