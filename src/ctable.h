#ifndef CAMIL_CTABLE_H
#define CAMIL_CTABLE_H


#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "logger.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// occupancy the table is sized for when it is frozen
#define CTABLE_TARGET_LOAD 0.55

// how many labels ahead the classification loop prefetches
#define CTABLE_PREFETCH_DISTANCE 8


struct ctable {
	struct camil_entry *entries; // capacity slots, one allocation
	uint64_t capacity;           // power of two
	uint64_t mask;               // capacity - 1
	uint64_t size;               // occupied slots
};

// allocates a table able to hold nentries entries at the target load factor
int ctable_init(struct ctable *table, uint64_t nentries);

// releases the table
void ctable_free(struct ctable *table);

// inserts entry, which must not already be present
int ctable_insert(struct ctable *table, const struct camil_entry *entry);

// number of entries stored, and the heap footprint in bytes
uint64_t ctable_size(const struct ctable *table);
uint64_t ctable_memory(const struct ctable *table);

// writes the raw entry array to out and reads it back
int ctable_read(struct ctable *table, FILE *in);
int ctable_write(const struct ctable *table, FILE *out);

// checks a table that came from a file
int ctable_validate(const struct ctable *table, uint32_t ngenomes);

// hints the memory system to start fetching the line that label would land on
static inline void ctable_prefetch(const struct ctable *table, lcp_label label) {
#if defined(__GNUC__)
	__builtin_prefetch(&table->entries[camil_hash(label) & table->mask], 0, 3);
#else
	(void)table;
	(void)label;
#endif
}

// returns the entry for label, or NULL when the label is absent
static inline const struct camil_entry *ctable_lookup(const struct ctable *table, lcp_label label) {
	uint64_t idx = camil_hash(label) & table->mask;

	for (;;) {
		const struct camil_entry *entry = &table->entries[idx];

		if (entry->ngenomes == 0) {
			return NULL; // empty slot ends the probe sequence
		}
		if (entry->label == label) {
			return entry;
		}
		idx = (idx + 1) & table->mask;
	}
}

#ifdef __cplusplus
}
#endif

#endif
