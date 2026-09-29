#include "cfull.h"
#include "cbuild.h"
#include "index.h"
// khashl compares an int loop bound with the unsigned ensemble width inside
// its own macros; nothing to fix on our side
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#include "khashl.h" // from deps


// The map is an ensemble of khashl tables keyed by the 64 bit core key. The
// value is the packed (offset, count) pair, see cfull.h. Everything khashl
// generates is file local; the header only forward declares struct cfull_map.
KHASHE_MAP_INIT(KH_LOCAL, cfull_map, cfull_map, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)

typedef cfull_map_sub cfull_sub;
typedef cfull_map_em_bucket_t cfull_bucket;

// bucket count of a khashl table, zero when it has never been sized
static inline uint64_t sub_capacity(const cfull_sub *sub) {
	return sub->keys != NULL ? (uint64_t)1 << sub->bits : 0;
}

// words in the used bitmap of a table with `capacity` buckets
static inline uint64_t sub_flag_words(uint64_t capacity) {
	return capacity < 32 ? 1 : capacity >> 5;
}

// number of sub tables the ensemble should be split into so that each holds
// about CFULL_SUB_TARGET entries
static int cfull_ensemble_bits(uint64_t ncores) {
	int bits = 0;

	while (bits < CFULL_MAX_ENSEMBLE_BITS && (ncores >> bits) > CFULL_SUB_TARGET) {
		bits++;
	}
	return bits;
}

int cfull_alloc(struct cfull *full, uint64_t ncores, uint64_t npairs) {
	int bits;
	uint64_t per_sub;
	uint64_t s;

	if (npairs > CFULL_MAX_PAIRS) {
		log_error("full index: %llu core occurrences exceed what the format can address", (unsigned long long)npairs);
		return -1;
	}

	bits = cfull_ensemble_bits(ncores);
	full->map = cfull_map_init(bits);
	full->sids = (camil_sid *)malloc((npairs > 0 ? npairs : 1) * sizeof(camil_sid));
	if (full->map == NULL || full->map->sub == NULL || full->sids == NULL) {
		log_error("full index: out of memory for %llu cores", (unsigned long long)ncores);
		cfull_map_destroy(full->map);
		free(full->sids);
		full->map = NULL;
		full->sids = NULL;
		return -1;
	}

	// Pre-size every sub table so that inserting never rehashes. The keys
	// spread evenly across the sub tables, being hashes, and a table that
	// still overflows simply grows on its own.
	per_sub = ncores >> bits;
	for (s = 0; s < ((uint64_t)1 << bits); s++) {
		uint64_t want = per_sub * 4 / 3 + 2;

		if (want > (uint64_t)1 << 31) {
			want = (uint64_t)1 << 31;
		}
		if (cfull_map_em_sub_resize(&full->map->sub[s], (khint_t)want) < 0) {
			log_error("full index: out of memory for %llu cores", (unsigned long long)ncores);
			cfull_map_destroy(full->map);
			free(full->sids);
			full->map = NULL;
			full->sids = NULL;
			return -1;
		}
	}

	full->ncores = ncores;
	full->nsids = 0;
	return 0;
}

int cfull_put(struct cfull *full, uint64_t key, uint64_t offset, uint32_t count) {
	kh_ensitr_t at;
	int absent;

	if (count > CFULL_MAX_COUNT) {
		log_error("full index: a core occurs in %u species, more than the format can record", count);
		return -1;
	}

	at = cfull_map_put(full->map, key, &absent);
	if (absent < 0) {
		log_error("full index: out of memory while inserting cores");
		return -1;
	}
	if (absent == 0) {
		log_error("full index: core %016llx was inserted twice", (unsigned long long)key);
		return -1;
	}
	kh_ens_val(full->map, at) = cfull_pack(offset, count);
	return 0;
}

int cfull_get(const struct cfull *full, uint64_t key, uint64_t *offset, uint32_t *count) {
	kh_ensitr_t at = cfull_map_get(full->map, key);
	uint64_t value;

	if (kh_ens_is_end(at)) {
		return -1;
	}
	value = kh_ens_val(full->map, at);
	*offset = cfull_offset(value);
	*count = cfull_count(value);
	return 0;
}

uint64_t cfull_size(const struct cfull *full) {
	return full->map != NULL ? (uint64_t)kh_ens_size(full->map) : 0;
}

uint64_t cfull_memory(const struct cfull *full) {
	uint64_t total = sizeof(struct cfull) + full->nsids * sizeof(camil_sid);
	uint64_t s;

	if (full->map != NULL) {
		for (s = 0; s < ((uint64_t)1 << full->map->bits); s++) {
			uint64_t capacity = sub_capacity(&full->map->sub[s]);

			total += capacity * sizeof(cfull_bucket) + sub_flag_words(capacity) * sizeof(khint32_t);
		}
	}
	return total;
}

void cfull_destroy(struct cfull *full) {
	uint32_t i;

	if (full->names != NULL) {
		for (i = 0; i < full->ngenomes; i++) {
			free(full->names[i]);
		}
		free(full->names);
	}
	cfull_map_destroy(full->map);
	free(full->sids);
	memset(full, 0, sizeof(*full));
}

const char *cfull_genome_name(const struct cfull *full, camil_sid sid) {
	if (full->names == NULL || sid >= full->ngenomes) {
		return "?";
	}
	return full->names[sid];
}

// -- serialization ----------------------------------------------------------

static int write_u32(FILE *out, uint32_t value) {
	return fwrite(&value, sizeof(value), 1, out) == 1 ? 0 : -1;
}

static int write_u64(FILE *out, uint64_t value) {
	return fwrite(&value, sizeof(value), 1, out) == 1 ? 0 : -1;
}

static int read_u32(FILE *in, uint32_t *value) {
	return fread(value, sizeof(*value), 1, in) == 1 ? 0 : -1;
}

static int read_u64(FILE *in, uint64_t *value) {
	return fread(value, sizeof(*value), 1, in) == 1 ? 0 : -1;
}

// dumps one khashl table: its shape, then its bitmap and buckets as they sit
// in memory
static int sub_write(const cfull_sub *sub, FILE *out) {
	uint64_t capacity = sub_capacity(sub);
	uint32_t present = capacity > 0;

	if (write_u32(out, sub->bits) != 0 || write_u32(out, sub->count) != 0 || write_u32(out, present) != 0) {
		return -1;
	}
	if (!present) {
		return 0;
	}
	if (fwrite(sub->used, sizeof(khint32_t), sub_flag_words(capacity), out) != sub_flag_words(capacity)) {
		return -1;
	}
	if (fwrite(sub->keys, sizeof(cfull_bucket), capacity, out) != capacity) {
		return -1;
	}
	return 0;
}

static int sub_read(cfull_sub *sub, FILE *in) {
	uint32_t bits;
	uint32_t count;
	uint32_t present;
	uint64_t capacity;

	if (read_u32(in, &bits) != 0 || read_u32(in, &count) != 0 || read_u32(in, &present) != 0) {
		return -1;
	}
	if (bits > 31 || present > 1) {
		return -1;
	}

	sub->bits = bits;
	sub->count = count;
	if (!present) {
		return count == 0 ? 0 : -1;
	}

	capacity = (uint64_t)1 << bits;
	if (count > capacity) {
		return -1;
	}
	sub->used = (khint32_t *)malloc(sub_flag_words(capacity) * sizeof(khint32_t));
	sub->keys = (cfull_bucket *)malloc(capacity * sizeof(cfull_bucket));
	if (sub->used == NULL || sub->keys == NULL) {
		log_error("full index: out of memory while loading");
		return -1;
	}
	if (fread(sub->used, sizeof(khint32_t), sub_flag_words(capacity), in) != sub_flag_words(capacity)) {
		return -1;
	}
	if (fread(sub->keys, sizeof(cfull_bucket), capacity, in) != capacity) {
		return -1;
	}
	return 0;
}

int cfull_save(const struct cfull *full, const char *path) {
	FILE *out;
	uint32_t i;
	uint64_t s;

	out = fopen(path, "wb");
	if (out == NULL) {
		log_error("cannot write the full index to %s", path);
		return -1;
	}

	if (fwrite(CFULL_MAGIC, 1, 8, out) != 8 || write_u32(out, CFULL_VERSION) != 0 || write_u32(out, 64) != 0) {
		goto write_error;
	}
	if (write_u32(out, (uint32_t)full->lcp_level) != 0 || write_u32(out, (uint32_t)full->use_rc) != 0 || write_u32(out, full->ngenomes) != 0) {
		goto write_error;
	}
	for (i = 0; i < full->ngenomes; i++) {
		uint32_t len = (uint32_t)strlen(full->names[i]);

		if (write_u32(out, len) != 0 || fwrite(full->names[i], 1, len, out) != len) {
			goto write_error;
		}
	}

	if (write_u64(out, cfull_size(full)) != 0 || write_u64(out, full->nsids) != 0) {
		goto write_error;
	}
	if (write_u32(out, (uint32_t)full->map->bits) != 0) {
		goto write_error;
	}
	for (s = 0; s < ((uint64_t)1 << full->map->bits); s++) {
		if (sub_write(&full->map->sub[s], out) != 0) {
			goto write_error;
		}
	}
	if (full->nsids > 0 && fwrite(full->sids, sizeof(camil_sid), full->nsids, out) != full->nsids) {
		goto write_error;
	}

	if (fclose(out) != 0) {
		log_error("cannot flush the full index to %s", path);
		return -1;
	}

	log_info("wrote %llu cores over %llu species occurrences to %s", (unsigned long long)cfull_size(full), (unsigned long long)full->nsids, path);
	return 0;

write_error:
	log_error("failed while writing the full index to %s", path);
	fclose(out);
	return -1;
}

int cfull_load(struct cfull *full, const char *path) {
	FILE *in;
	char magic[8];
	uint32_t value;
	uint32_t g;
	uint64_t ncores;
	uint64_t counted = 0;
	uint64_t s;

	memset(full, 0, sizeof(*full));

	in = fopen(path, "rb");
	if (in == NULL) {
		log_error("cannot open the full index %s", path);
		return -1;
	}

	if (fread(magic, 1, 8, in) != 8 || memcmp(magic, CFULL_MAGIC, 8) != 0) {
		log_error("%s is not a camil full index", path);
		goto read_error;
	}
	if (read_u32(in, &value) != 0 || value != CFULL_VERSION) {
		log_error("%s was written in format version %u, this build reads version %u; rebuild the index", path, value, CFULL_VERSION);
		goto read_error;
	}
	if (read_u32(in, &value) != 0 || value != 64) {
		log_error("%s stores %u bit keys but this build uses 64 bit keys", path, value);
		goto read_error;
	}

	if (read_u32(in, &value) != 0 || value < 1 || value > INT_MAX) {
		goto truncated;
	}
	full->lcp_level = (int)value;
	if (read_u32(in, &value) != 0) {
		goto truncated;
	}
	full->use_rc = (int)value;
	if (read_u32(in, &value) != 0 || value == 0 || value > CAMIL_MAX_GENOMES) {
		goto truncated;
	}
	full->ngenomes = value;

	full->names = (char **)calloc(full->ngenomes, sizeof(char *));
	if (full->names == NULL) {
		log_error("out of memory while loading %s", path);
		goto read_error;
	}
	for (g = 0; g < full->ngenomes; g++) {
		uint32_t len;

		if (read_u32(in, &len) != 0 || len > (1u << 16)) {
			goto truncated;
		}
		full->names[g] = (char *)malloc(len + 1);
		if (full->names[g] == NULL) {
			log_error("out of memory while loading %s", path);
			goto read_error;
		}
		if (fread(full->names[g], 1, len, in) != len) {
			goto truncated;
		}
		full->names[g][len] = '\0';
	}

	if (read_u64(in, &ncores) != 0 || read_u64(in, &full->nsids) != 0 || full->nsids > CFULL_MAX_PAIRS) {
		goto truncated;
	}
	if (read_u32(in, &value) != 0 || value > CFULL_MAX_ENSEMBLE_BITS) {
		goto truncated;
	}

	full->map = cfull_map_init((int)value);
	if (full->map == NULL || full->map->sub == NULL) {
		log_error("out of memory while loading %s", path);
		goto read_error;
	}
	for (s = 0; s < ((uint64_t)1 << full->map->bits); s++) {
		if (sub_read(&full->map->sub[s], in) != 0) {
			goto truncated;
		}
		counted += full->map->sub[s].count;
	}
	if (counted != ncores) {
		log_error("%s is corrupted: %llu cores announced, %llu stored", path, (unsigned long long)ncores, (unsigned long long)counted);
		goto read_error;
	}
	full->map->count = counted;
	full->ncores = ncores;

	full->sids = (camil_sid *)malloc((full->nsids > 0 ? full->nsids : 1) * sizeof(camil_sid));
	if (full->sids == NULL) {
		log_error("out of memory while loading %s", path);
		goto read_error;
	}
	if (full->nsids > 0 && fread(full->sids, sizeof(camil_sid), full->nsids, in) != full->nsids) {
		goto truncated;
	}

	// Everything above came from a file, and the subset extraction indexes
	// the name table and the species array with what it finds in the map.
	for (s = 0; s < full->nsids; s++) {
		if (full->sids[s] >= full->ngenomes) {
			log_error("%s is corrupted: an entry refers to species %u of %u", path, full->sids[s], full->ngenomes);
			goto read_error;
		}
	}
	{
		kh_ensitr_t it;

		kh_ens_foreach(full->map, it) {
			uint64_t v = kh_ens_val(full->map, it);

			if (cfull_count(v) == 0 || cfull_offset(v) + cfull_count(v) > full->nsids) {
				log_error("%s is corrupted: a core points outside the species array", path);
				goto read_error;
			}
		}
	}

	fclose(in);
	log_info("loaded %llu cores over %llu species occurrences from %s", (unsigned long long)ncores, (unsigned long long)full->nsids, path);
	return 0;

truncated:
	log_error("%s is truncated or corrupted", path);

read_error:
	fclose(in);
	cfull_destroy(full);
	return -1;
}

void cfull_report(const struct cfull *full) {
	uint32_t i;

	log_info("full index: LCP level %d, reverse complement %s, %u species", full->lcp_level, full->use_rc ? "on" : "off", full->ngenomes);
	if (full->ngenomes <= 64) {
		for (i = 0; i < full->ngenomes; i++) {
			log_info("  species %u: %s", i, full->names[i]);
		}
	}
	log_info("full index: %llu cores, %llu species occurrences, about %llu MiB in memory", (unsigned long long)cfull_size(full), (unsigned long long)full->nsids, (unsigned long long)(cfull_memory(full) / (1024ull * 1024ull)));
}

// -- subset extraction ------------------------------------------------------

// species id in `full` that carries `name`, or -1 when none does
static int64_t cfull_find_name(const struct cfull *full, const char *name) {
	uint32_t i;

	for (i = 0; i < full->ngenomes; i++) {
		if (strcmp(full->names[i], name) == 0) {
			return (int64_t)i;
		}
	}
	return -1;
}

int cfull_subset(const struct cfull *full, char *const *names, uint32_t nnames, uint32_t max_share, int threads, struct camil_index *index) {
	camil_sid *renumber = NULL; // old species id -> new, or CAMIL_SID_NONE
	struct cbuild builder;
	struct cbuild_sink *sink;
	kh_ensitr_t it;
	uint64_t kept = 0;
	uint64_t dropped_shared = 0;
	uint32_t i;
	int status = 0;

	memset(index, 0, sizeof(*index));
	memset(&builder, 0, sizeof(builder));

	if (nnames == 0) {
		log_error("subset: at least one species name is required");
		return -1;
	}
	if (max_share < 1) {
		log_error("the sharing limit must be at least 1");
		return -1;
	}

	renumber = (camil_sid *)malloc((full->ngenomes > 0 ? full->ngenomes : 1) * sizeof(camil_sid));
	index->names = (char **)calloc(nnames, sizeof(char *));
	if (renumber == NULL || index->names == NULL) {
		log_error("out of memory while preparing the subset");
		status = -1;
		goto done;
	}
	for (i = 0; i < full->ngenomes; i++) {
		renumber[i] = CAMIL_SID_NONE;
	}

	for (i = 0; i < nnames; i++) {
		int64_t old = cfull_find_name(full, names[i]);

		if (old < 0) {
			log_error("subset: the index has no species named '%s'", names[i]);
			status = -1;
			goto done;
		}
		if (renumber[old] != CAMIL_SID_NONE) {
			log_error("subset: species '%s' was named twice", names[i]);
			status = -1;
			goto done;
		}
		renumber[old] = (camil_sid)i;
		index->names[i] = camil_strdup(full->names[old]);
		if (index->names[i] == NULL) {
			log_error("out of memory while preparing the subset");
			status = -1;
			goto done;
		}
		log_info("species %u: %s (was %lld)", i, index->names[i], (long long)old);
	}
	index->ngenomes = nnames;
	index->lcp_level = full->lcp_level;
	index->max_share = (int)max_share;
	index->use_rc = full->use_rc;

	// The kept pairs go through the same sort as a fresh build, which is what
	// turns the hash order of the map into the sorted table classify wants.
	// One sink is enough: the walk over the map is sequential and the sort
	// behind it is what runs in parallel.
	if (cbuild_init(&builder, threads) != 0) {
		status = -1;
		goto done;
	}
	sink = &builder.sinks[0];

	kh_ens_foreach(full->map, it) {
		uint64_t value = kh_ens_val(full->map, it);
		uint64_t key = kh_ens_key(full->map, it);
		const camil_sid *run = full->sids + cfull_offset(value);
		uint32_t count = cfull_count(value);
		uint32_t chosen = 0;
		uint32_t k;

		for (k = 0; k < count; k++) {
			chosen += renumber[run[k]] != CAMIL_SID_NONE;
		}
		if (chosen == 0) {
			continue;
		}
		if (chosen > max_share) {
			dropped_shared++;
			continue;
		}
		for (k = 0; k < count; k++) {
			camil_sid sid = renumber[run[k]];

			if (sid != CAMIL_SID_NONE && cbuild_push(&builder, sink, key, sid) != 0) {
				status = -1;
				goto done;
			}
		}
		kept++;
	}

	log_info("subset: %llu cores occur in the %u chosen species, %llu of them in more than %u and dropped", (unsigned long long)(kept + dropped_shared), nnames, (unsigned long long)dropped_shared, max_share);

	// nothing is filtered here: the sharing limit was applied above
	if (cbuild_freeze(&builder, CBUILD_KEEP_ALL, threads, &index->table) != 0) {
		status = -1;
		goto done;
	}

done:
	cbuild_free(&builder);
	free(renumber);
	if (status != 0) {
		camil_index_destroy(index);
	}
	return status;
}

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
