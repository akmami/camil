#include "ctable.h"


// smallest table, so that an empty index is still a valid one to query
#define CTABLE_MIN_CAPACITY 16u

static uint64_t ctable_next_pow2(uint64_t value) {
	uint64_t result = CTABLE_MIN_CAPACITY;

	while (result < value) {
		result <<= 1;
	}
	return result;
}

int ctable_init(struct ctable *table, uint64_t nentries) {
	uint64_t capacity = ctable_next_pow2((uint64_t)((double)nentries / CTABLE_TARGET_LOAD) + 1);

	memset(table, 0, sizeof(*table));

	table->entries = (struct camil_entry *)calloc(capacity, sizeof(struct camil_entry));
	if (table->entries == NULL) {
		log_error("core table: out of memory while allocating %llu slots (%llu MiB)", (unsigned long long)capacity, (unsigned long long)(capacity * sizeof(struct camil_entry) / (1024 * 1024)));
		return -1;
	}
	table->capacity = capacity;
	table->mask = capacity - 1;
	table->size = 0;
	return 0;
}

void ctable_free(struct ctable *table) {
	free(table->entries);
	memset(table, 0, sizeof(*table));
}

int ctable_insert(struct ctable *table, const struct camil_entry *entry) {
	uint64_t idx = camil_hash(entry->label) & table->mask;

	if (table->size + 1 >= table->capacity) {
		log_error("core table: no room left for %llu entries in %llu slots", (unsigned long long)table->size + 1, (unsigned long long)table->capacity);
		return -1;
	}

	while (table->entries[idx].ngenomes != 0) {
		idx = (idx + 1) & table->mask;
	}
	table->entries[idx] = *entry;
	table->size++;
	return 0;
}

uint64_t ctable_size(const struct ctable *table) {
	return table->size;
}

uint64_t ctable_memory(const struct ctable *table) {
	return table->capacity * sizeof(struct camil_entry) + sizeof(struct ctable);
}

int ctable_write(const struct ctable *table, FILE *out) {
	if (fwrite(&table->capacity, sizeof(table->capacity), 1, out) != 1 ||
	    fwrite(&table->size, sizeof(table->size), 1, out) != 1) {
		return -1;
	}
	if (table->capacity > 0 &&
	    fwrite(table->entries, sizeof(struct camil_entry), table->capacity, out) != table->capacity) {
		return -1;
	}
	return 0;
}

int ctable_read(struct ctable *table, FILE *in) {
	uint64_t capacity;
	uint64_t size;

	memset(table, 0, sizeof(*table));

	if (fread(&capacity, sizeof(capacity), 1, in) != 1 ||
	    fread(&size, sizeof(size), 1, in) != 1) {
		return -1;
	}

	// capacity is a power of two by construction
	if (capacity == 0 || (capacity & (capacity - 1)) != 0 || size >= capacity) {
		log_error("core table: the stored table has an invalid shape");
		return -1;
	}

	// capacity comes straight from the file
	if (capacity > SIZE_MAX / sizeof(struct camil_entry)) {
		log_error("core table: the stored table claims an impossible size");
		return -1;
	}

	table->entries = (struct camil_entry *)malloc(capacity * sizeof(struct camil_entry));
	if (table->entries == NULL) {
		log_error("core table: out of memory while loading %llu slots", (unsigned long long)capacity);
		return -1;
	}
	if (fread(table->entries, sizeof(struct camil_entry), capacity, in) != capacity) {
		free(table->entries);
		memset(table, 0, sizeof(*table));
		return -1;
	}

	table->capacity = capacity;
	table->mask = capacity - 1;
	table->size = size;
	return 0;
}

int ctable_validate(const struct ctable *table, uint32_t ngenomes) {
	uint64_t occupied = 0;
	uint64_t i;

	// everything in the slot array arrives unchecked from the file
	// classification loop indexes a stack array with the genome ids it finds here
	for (i = 0; i < table->capacity; i++) {
		const struct camil_entry *entry = &table->entries[i];
		uint8_t k;

		if (entry->ngenomes == 0) {
			continue;
		}
		if (entry->ngenomes > CAMIL_MAX_SHARE) {
			log_error("core table: an entry claims %u genomes, the limit is %d", entry->ngenomes, CAMIL_MAX_SHARE);
			return -1;
		}
		for (k = 0; k < entry->ngenomes; k++) {
			if (entry->gids[k] >= ngenomes) {
				log_error("core table: an entry refers to genome %u of %u", entry->gids[k], ngenomes);
				return -1;
			}
		}
		occupied++;
	}

	if (occupied != table->size) {
		log_error("core table: %llu slots are occupied but the header says %llu", (unsigned long long)occupied, (unsigned long long)table->size);
		return -1;
	}

	return 0;
}
