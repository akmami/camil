#include "classify.h"


// outcome for a single read
struct read_result {
	uint32_t ncores;   // cores the read decomposed into
	uint32_t nmatched; // cores found in the index
	uint16_t best;     // votes for the winning genome, saturating
	uint16_t second;   // votes for the runner up, saturating
	int16_t gid;       // winning genome, or -1 when there is none
	uint8_t status;    // enum camil_status
	uint8_t reserved;  // keeps the record at 16 bytes
};

// number of results that share one cache line; chunk sizes are rounded to it
#define RESULTS_PER_LINE (CAMIL_CACHE_LINE / (int)sizeof(struct read_result))

// batch of reads plus the space their results go into
struct read_batch {
	char **names;                // pointers into the name arena
	char **seqs;                 // pointers into the sequence arena
	uint32_t *lens;              // sequence lengths
	struct read_result *results; // one per read
	uint32_t nreads;             // reads currently held
	uint32_t capacity;           // reads the arrays can hold

	char *name_arena; // concatenated NUL terminated names
	size_t name_used;
	size_t name_capacity;
	char *seq_arena;  // concatenated NUL terminated sequences
	size_t seq_used;
	size_t seq_capacity;
};

// one chunk of a batch handed to a worker
struct classify_job {
	const struct camil_index *index;
	const struct camil_thresholds *thresholds;
	struct read_batch *batch;
	uint32_t begin; // first read of the chunk
	uint32_t end;   // one past the last read of the chunk
};

void camil_thresholds_default(struct camil_thresholds *thresholds) {
	thresholds->min_hits = 1;
	thresholds->min_ratio = 1.0;
}

const char *camil_status_name(int status) {
	switch (status) {
	case CAMIL_ASSIGNED:
		return "assigned";
	case CAMIL_AMBIGUOUS:
		return "ambiguous";
	case CAMIL_UNCLASSIFIED:
		return "unclassified";
	default:
		return "unknown";
	}
}

int camil_stats_init(struct camil_stats *stats, uint32_t ngenomes) {
	memset(stats, 0, sizeof(*stats));
	stats->assigned = (uint64_t *)calloc(ngenomes > 0 ? ngenomes : 1, sizeof(uint64_t));
	if (stats->assigned == NULL) {
		log_error("out of memory while allocating the report counters");
		return -1;
	}
	stats->ngenomes = ngenomes;
	return 0;
}

void camil_stats_destroy(struct camil_stats *stats) {
	free(stats->assigned);
	memset(stats, 0, sizeof(*stats));
}

void camil_stats_merge(struct camil_stats *dst, const struct camil_stats *src) {
	uint32_t i;

	for (i = 0; i < dst->ngenomes && i < src->ngenomes; i++) {
		dst->assigned[i] += src->assigned[i];
	}
	dst->ambiguous += src->ambiguous;
	dst->unclassified += src->unclassified;
	dst->total += src->total;
}

// grows an arena so that needed more bytes fit. Returns 0 on success.
static int arena_reserve(char **arena, size_t used, size_t *capacity, size_t needed) {
	size_t target = used + needed;
	size_t grown_capacity = *capacity;
	char *grown;

	if (target <= *capacity) {
		return 0;
	}
	while (grown_capacity < target) {
		grown_capacity = grown_capacity == 0 ? (1u << 16) : grown_capacity * 2;
	}
	// new capacity is only recorded once the memory behind it exists,
	// otherwise a failed reallocation would leave the arena describing space
	// it does not own.
	grown = (char *)realloc(*arena, grown_capacity);
	if (grown == NULL) {
		log_error("out of memory while buffering reads");
		return -1;
	}
	*arena = grown;
	*capacity = grown_capacity;
	return 0;
}

static int batch_init(struct read_batch *batch, uint32_t capacity) {
	memset(batch, 0, sizeof(*batch));
	batch->names = (char **)malloc(capacity * sizeof(char *));
	batch->seqs = (char **)malloc(capacity * sizeof(char *));
	batch->lens = (uint32_t *)malloc(capacity * sizeof(uint32_t));
	batch->results = (struct read_result *)malloc(capacity * sizeof(struct read_result));
	if (batch->names == NULL || batch->seqs == NULL || batch->lens == NULL ||
	    batch->results == NULL) {
		log_error("out of memory while allocating a read batch");
		return -1;
	}
	batch->capacity = capacity;
	return 0;
}

static void batch_destroy(struct read_batch *batch) {
	free(batch->names);
	free(batch->seqs);
	free(batch->lens);
	free(batch->results);
	free(batch->name_arena);
	free(batch->seq_arena);
	memset(batch, 0, sizeof(*batch));
}

static void batch_clear(struct read_batch *batch) {
	batch->nreads = 0;
	batch->name_used = 0;
	batch->seq_used = 0;
}

// copies one record into the batch
static int batch_push(struct read_batch *batch, const struct seqrec *rec) {
	size_t name_len = strlen(rec->name);
	size_t name_offset = batch->name_used;
	size_t seq_offset = batch->seq_used;
	char *old_names = batch->name_arena;
	char *old_seqs = batch->seq_arena;
	int status = 0;
	uint32_t i;

	if (arena_reserve(&batch->name_arena, batch->name_used, &batch->name_capacity, name_len + 1) != 0) {
		status = -1;
	}
	if (status == 0 && arena_reserve(&batch->seq_arena, batch->seq_used, &batch->seq_capacity, rec->len + 1) != 0) {
		status = -1;
	}

	// rebase the pointers of the reads already in the batch if an arena moved
	if (batch->name_arena != old_names || batch->seq_arena != old_seqs) {
		for (i = 0; i < batch->nreads; i++) {
			batch->names[i] = batch->name_arena + (size_t)(batch->names[i] - old_names);
			batch->seqs[i] = batch->seq_arena + (size_t)(batch->seqs[i] - old_seqs);
		}
	}

	if (status != 0) {
		return -1;
	}

	memcpy(batch->name_arena + name_offset, rec->name, name_len + 1);
	memcpy(batch->seq_arena + seq_offset, rec->seq, rec->len + 1);
	batch->name_used += name_len + 1;
	batch->seq_used += rec->len + 1;

	batch->names[batch->nreads] = batch->name_arena + name_offset;
	batch->seqs[batch->nreads] = batch->seq_arena + seq_offset;
	batch->lens[batch->nreads] = (uint32_t)rec->len;
	batch->nreads++;
	return 0;
}

// clamps a vote count to what the report record can hold
static uint16_t saturate16(uint32_t value) {
	return value > 0xFFFFu ? 0xFFFFu : (uint16_t)value;
}

// scores one read: parse, look every core up, tally votes, apply thresholds
static void classify_read(const struct camil_index *index, const struct camil_thresholds *thresholds, const char *seq, uint32_t len, uint32_t *votes, struct read_result *result) {
	const struct ctable *table = &index->table;
	struct lps parsed;
	uint32_t nmatched = 0;
	uint32_t best = 0;
	uint32_t second = 0;
	int32_t best_gid = -1;
	int tied = 0;
	int i;
	int n;
	uint32_t g;

	memset(result, 0, sizeof(*result));
	result->gid = -1;
	result->status = CAMIL_UNCLASSIFIED;

	if (len == 0) {
		return;
	}

	memset(votes, 0, index->ngenomes * sizeof(uint32_t));

	init_lps(&parsed, seq, (int)len);
	lps_deepen(&parsed, index->lcp_level);
	n = parsed.size;

	// warm up the pipeline: the first few slots are requested before the first lookup so that the loop below always has a fetch in flight
	for (i = 0; i < CTABLE_PREFETCH_DISTANCE && i < n; i++) {
		ctable_prefetch(table, parsed.cores[i].label);
	}

	for (i = 0; i < n; i++) {
		const struct camil_entry *entry;
		uint8_t k;

		if (i + CTABLE_PREFETCH_DISTANCE < n) {
			ctable_prefetch(table, parsed.cores[i + CTABLE_PREFETCH_DISTANCE].label);
		}

		entry = ctable_lookup(table, parsed.cores[i].label);
		if (entry == NULL) {
			continue;
		}
		nmatched++;
		for (k = 0; k < entry->ngenomes; k++) {
			votes[entry->gids[k]]++;
		}
	}

	result->ncores = n > 0 ? (uint32_t)n : 0;
	free_lps(&parsed);

	result->nmatched = nmatched;
	if (nmatched == 0) {
		return;
	}

	// highest and second highest vote count, remembering whether the top is shared by more than one genome
	for (g = 0; g < index->ngenomes; g++) {
		if (votes[g] > best) {
			second = best;
			best = votes[g];
			best_gid = (int32_t)g;
			tied = 0;
		} else if (votes[g] == best && best > 0) {
			tied = 1;
			second = best;
		} else if (votes[g] > second) {
			second = votes[g];
		}
	}

	result->best = saturate16(best);
	result->second = saturate16(second);

	if (tied || best == 0) {
		result->status = CAMIL_AMBIGUOUS;
		return;
	}
	if (best < thresholds->min_hits) {
		result->status = CAMIL_AMBIGUOUS;
		return;
	}
	if ((double)best < thresholds->min_ratio * (double)nmatched) {
		result->status = CAMIL_AMBIGUOUS;
		return;
	}

	result->gid = (int16_t)best_gid;
	result->status = CAMIL_ASSIGNED;
}

// thread pool entry point: classify one contiguous chunk of a batch.
static void classify_job_run(void *arg) {
	struct classify_job *job = (struct classify_job *)arg;
	// Thread private vote counters
	uint32_t votes[CAMIL_MAX_GENOMES];
	uint32_t i;

	for (i = job->begin; i < job->end; i++) {
		classify_read(job->index, job->thresholds, job->batch->seqs[i], job->batch->lens[i], votes, &job->batch->results[i]);
	}
}

void camil_classify_write_header(FILE *per_read) {
	if (per_read == NULL) {
		return;
	}
	fprintf(per_read, "#read\tlength\tcores\tmatched\tbest_hits\tsecond_hits\tassignment\tstatus\n");
}

// writes the per read lines of a finished batch and folds them into the tally
static void batch_report(const struct read_batch *batch, const struct camil_index *index, FILE *per_read, struct camil_stats *stats) {
	uint32_t i;

	for (i = 0; i < batch->nreads; i++) {
		const struct read_result *result = &batch->results[i];
		const char *assignment = "-";

		if (result->status == CAMIL_ASSIGNED) {
			assignment = camil_index_genome_name(index, (camil_gid)result->gid);
			stats->assigned[result->gid]++;
		} else if (result->status == CAMIL_AMBIGUOUS) {
			stats->ambiguous++;
		} else {
			stats->unclassified++;
		}
		stats->total++;

		if (per_read != NULL) {
			fprintf(per_read, "%s\t%u\t%u\t%u\t%u\t%u\t%s\t%s\n",
			        batch->names[i], batch->lens[i], result->ncores, result->nmatched,
			        result->best, result->second, assignment,
			        camil_status_name(result->status));
		}
	}
}

// fills batch with the next records. Returns 0 on success and -1 on a malformed file or an allocation failure; an empty batch means end of input
static int batch_fill(struct read_batch *batch, struct seqfile *file, const char *path) {
	struct seqrec rec;
	int ret = 0;

	batch_clear(batch);

	while (batch->nreads < batch->capacity && (ret = seq_read(file, &rec)) > 0) {
		// lcptools parses a sequence in one piece; a read longer than the
		// parser limit is reported rather than silently mishandled.
		if (rec.len > (size_t)INT_MAX) {
			log_warn("%s: read %s is too long to parse, skipped", path, rec.name);
			continue;
		}
		if (batch_push(batch, &rec) != 0) {
			return -1;
		}
	}

	return ret < 0 ? -1 : 0;
}

// splits batch into chunks and hands them to the pool
static void batch_dispatch(struct read_batch *batch, struct classify_job *jobs, int nchunks, const struct camil_index *index, const struct camil_thresholds *thresholds, struct tpool *pool) {
	uint32_t per_chunk;
	uint32_t cursor = 0;
	int j;

	per_chunk = (batch->nreads + (uint32_t)nchunks - 1) / (uint32_t)nchunks;
	per_chunk = ((per_chunk + RESULTS_PER_LINE - 1) / RESULTS_PER_LINE) * RESULTS_PER_LINE;
	if (per_chunk == 0) {
		per_chunk = RESULTS_PER_LINE;
	}

	for (j = 0; j < nchunks && cursor < batch->nreads; j++) {
		jobs[j].index = index;
		jobs[j].thresholds = thresholds;
		jobs[j].batch = batch;
		jobs[j].begin = cursor;
		jobs[j].end = cursor + per_chunk;
		if (jobs[j].end > batch->nreads) {
			jobs[j].end = batch->nreads;
		}
		cursor = jobs[j].end;
		tpool_submit(pool, classify_job_run, &jobs[j]);
	}
}

int camil_classify_file(const struct camil_index *index, const char *path, const struct camil_thresholds *thresholds, int threads, FILE *per_read, struct camil_stats *stats) {
	struct seqfile *file;
	struct read_batch batches[2];
	struct tpool pool;
	struct classify_job *jobs = NULL;
	int nchunks = threads > 1 ? threads : 1;
	int current = 0;
	int status = 0;
	int i;

	file = seq_open(path);
	if (file == NULL) {
		return -1;
	}

	// two batches: one being classified, one being filled
	for (i = 0; i < 2; i++) {
		if (batch_init(&batches[i], CAMIL_BATCH_READS) != 0) {
			status = -1;
		}
	}

	jobs = (struct classify_job *)calloc((size_t)nchunks, sizeof(struct classify_job));
	if (jobs == NULL) {
		log_error("out of memory while preparing the classification workers");
		status = -1;
	}

	if (status == 0 && tpool_init(&pool, threads, (size_t)nchunks) != 0) {
		status = -1;
	}

	if (status == 0) {
		status = batch_fill(&batches[current], file, path);

		while (status == 0 && batches[current].nreads > 0) {
			int next = current ^ 1;

			batch_dispatch(&batches[current], jobs, nchunks, index, thresholds, &pool);

			// read the next batch while the workers are busy with this one
			status = batch_fill(&batches[next], file, path);

			tpool_wait(&pool);
			batch_report(&batches[current], index, per_read, stats);
			current = next;
		}

		tpool_destroy(&pool);
	}

	free(jobs);
	for (i = 0; i < 2; i++) {
		batch_destroy(&batches[i]);
	}
	seq_close(file);

	return status;
}

void camil_stats_report(const struct camil_stats *stats, const struct camil_index *index, const char *label, FILE *out) {
	uint32_t i;
	double total = stats->total > 0 ? (double)stats->total : 1.0;

	fprintf(out, "# %s\n", label);
	fprintf(out, "#category\treads\tpercent\n");
	for (i = 0; i < stats->ngenomes; i++) {
		fprintf(out, "%s\t%llu\t%.2f\n", camil_index_genome_name(index, (camil_gid)i), (unsigned long long)stats->assigned[i], 100.0 * (double)stats->assigned[i] / total);
	}
	fprintf(out, "ambiguous\t%llu\t%.2f\n", (unsigned long long)stats->ambiguous, 100.0 * (double)stats->ambiguous / total);
	fprintf(out, "unclassified\t%llu\t%.2f\n", (unsigned long long)stats->unclassified, 100.0 * (double)stats->unclassified / total);
	fprintf(out, "total\t%llu\t100.00\n", (unsigned long long)stats->total);
	fflush(out);
}
