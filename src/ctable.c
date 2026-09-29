#include "ctable.h"


// Picks the directory size so that a bucket spans CTABLE_BUCKET_TARGET entries
// on average. Keys are hash values, so occupancy is even and the average is
// also close to the typical case.
static uint32_t ctable_dir_bits(uint64_t nentries) {
	uint64_t buckets = nentries / CTABLE_BUCKET_TARGET + 1;
	uint32_t bits = CTABLE_MIN_DIR_BITS;

	while (bits < CTABLE_MAX_DIR_BITS && ((uint64_t)1 << bits) < buckets) {
		bits++;
	}
	return bits;
}

int ctable_alloc(struct ctable *table, uint64_t nentries) {
	uint64_t ndir;

	memset(table, 0, sizeof(*table));

	// Directory offsets are 32 bit, which is what caps an index here. Four
	// billion cores is far past anything the rest of the tool handles, but the
	// check keeps the failure explicit rather than silently truncating.
	if (nentries > UINT32_MAX) {
		log_error("core table: %llu cores exceeds what one index can address", (unsigned long long)nentries);
		return -1;
	}

	table->dir_bits = ctable_dir_bits(nentries);
	ndir = ((uint64_t)1 << table->dir_bits) + 1;

	table->directory = (uint32_t *)calloc(ndir, sizeof(uint32_t));
	table->keys = (uint64_t *)malloc((nentries > 0 ? nentries : 1) * sizeof(uint64_t));
	table->values = (struct camil_value *)malloc((nentries > 0 ? nentries : 1) * sizeof(struct camil_value));

	if (table->directory == NULL || table->keys == NULL || table->values == NULL) {
		log_error("core table: out of memory for %llu cores (%llu MiB)", (unsigned long long)nentries, (unsigned long long)(nentries * 12 / (1024 * 1024)));
		ctable_free(table);
		return -1;
	}

	table->nentries = nentries;
	return 0;
}

void ctable_index(struct ctable *table) {
	uint64_t ndir = ((uint64_t)1 << table->dir_bits) + 1;
	uint64_t i;
	uint64_t bucket;
	uint64_t next = 0;

	// One pass writes, for every bucket, where its run starts. Buckets with no
	// keys end up equal to their neighbour, which makes their run empty.
	for (i = 0; i < table->nentries; i++) {
		bucket = ctable_bucket(table, table->keys[i]);
		while (next <= bucket) {
			table->directory[next++] = (uint32_t)i;
		}
	}
	while (next < ndir) {
		table->directory[next++] = (uint32_t)table->nentries;
	}
}

void ctable_free(struct ctable *table) {
	free(table->directory);
	free(table->keys);
	free(table->values);
	memset(table, 0, sizeof(*table));
}

uint64_t ctable_memory(const struct ctable *table) {
	uint64_t ndir = ((uint64_t)1 << table->dir_bits) + 1;

	return ndir * sizeof(uint32_t) + table->nentries * (sizeof(uint64_t) + sizeof(struct camil_value)) + sizeof(struct ctable);
}

int ctable_write(const struct ctable *table, FILE *out) {
	uint64_t ndir = ((uint64_t)1 << table->dir_bits) + 1;
	uint32_t bits = table->dir_bits;

	if (fwrite(&bits, sizeof(bits), 1, out) != 1 || fwrite(&table->nentries, sizeof(table->nentries), 1, out) != 1) {
		return -1;
	}
	if (fwrite(table->directory, sizeof(uint32_t), ndir, out) != ndir) {
		return -1;
	}
	if (table->nentries > 0) {
		if (fwrite(table->keys, sizeof(uint64_t), table->nentries, out) != table->nentries) {
			return -1;
		}
		if (fwrite(table->values, sizeof(struct camil_value), table->nentries, out) != table->nentries) {
			return -1;
		}
	}
	return 0;
}

int ctable_read(struct ctable *table, FILE *in) {
	uint32_t bits;
	uint64_t nentries;
	uint64_t ndir;

	memset(table, 0, sizeof(*table));

	if (fread(&bits, sizeof(bits), 1, in) != 1 || fread(&nentries, sizeof(nentries), 1, in) != 1) {
		return -1;
	}
	if (bits < CTABLE_MIN_DIR_BITS || bits > CTABLE_MAX_DIR_BITS || nentries > UINT32_MAX) {
		log_error("core table: the stored table has an invalid shape");
		return -1;
	}

	ndir = ((uint64_t)1 << bits) + 1;

	table->directory = (uint32_t *)malloc(ndir * sizeof(uint32_t));
	table->keys = (uint64_t *)malloc((nentries > 0 ? nentries : 1) * sizeof(uint64_t));
	table->values = (struct camil_value *)malloc((nentries > 0 ? nentries : 1) * sizeof(struct camil_value));
	if (table->directory == NULL || table->keys == NULL || table->values == NULL) {
		log_error("core table: out of memory while loading %llu cores", (unsigned long long)nentries);
		ctable_free(table);
		return -1;
	}

	table->dir_bits = bits;
	table->nentries = nentries;

	if (fread(table->directory, sizeof(uint32_t), ndir, in) != ndir) {
		ctable_free(table);
		return -1;
	}
	if (nentries > 0) {
		if (fread(table->keys, sizeof(uint64_t), nentries, in) != nentries || fread(table->values, sizeof(struct camil_value), nentries, in) != nentries) {
			ctable_free(table);
			return -1;
		}
	}
	return 0;
}

int ctable_validate(const struct ctable *table, uint32_t ngenomes) {
	uint64_t ndir = ((uint64_t)1 << table->dir_bits) + 1;
	uint32_t bucket;
	uint64_t i;

	// Everything below arrives unchecked from a file, and the classification
	// loop indexes a vote array with the species ids it finds here.
	if (table->directory[0] != 0 || table->directory[ndir - 1] != table->nentries) {
		log_error("core table: the directory does not span the entries");
		return -1;
	}
	for (i = 1; i < ndir; i++) {
		if (table->directory[i] < table->directory[i - 1] || table->directory[i] > table->nentries) {
			log_error("core table: the directory is not monotonic");
			return -1;
		}
	}
	for (i = 0; i < table->nentries; i++) {
		if (i > 0 && table->keys[i] < table->keys[i - 1]) {
			log_error("core table: the keys are not sorted");
			return -1;
		}
		bucket = ctable_bucket(table, table->keys[i]);
		if (i < table->directory[bucket] || i >= table->directory[bucket + 1]) {
			log_error("core table: the directory does not agree with the keys");
			return -1;
		}
		if (table->values[i].sid >= ngenomes) {
			log_error("core table: an entry refers to species %u of %u", table->values[i].sid, ngenomes);
			return -1;
		}
	}
	return 0;
}
