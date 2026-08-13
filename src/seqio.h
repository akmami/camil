#ifndef CAMIL_SEQIO_H
#define CAMIL_SEQIO_H


#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include "logger.h"
#include "kseq.h" // from deps
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stddef.h>


// opaque handle to an open sequence file
struct seqfile;

// single record
struct seqrec {
	const char *name;
	const char *seq;
	size_t len;
};

// opens path for reading, transparently handling gzip compressed files
struct seqfile *seq_open(const char *path);

// reads the next record into rec
int seq_read(struct seqfile *file, struct seqrec *rec);

// closes the file and releases the reader
void seq_close(struct seqfile *file);

// returns 1 when path exists and is readable, 0 otherwise
int seq_readable(const char *path);

// derives a display name from a file path
char *seq_basename(const char *path);

#ifdef __cplusplus
}
#endif

#endif
