#include "seqio.h"


KSEQ_INIT(gzFile, gzread)

struct seqfile {
	gzFile fp;      // compressed or plain input stream
	kseq_t *kseq;   // kseq parser state, owns the record buffers
	char *path;     // copy of the path, kept for error messages
};

struct seqfile *seq_open(const char *path) {
	struct seqfile *file;

	file = (struct seqfile *)calloc(1, sizeof(struct seqfile));
	if (file == NULL) {
		log_error("out of memory while opening %s", path);
		return NULL;
	}

	file->fp = gzopen(path, "r");
	if (file->fp == NULL) {
		log_error("cannot open %s", path);
		free(file);
		return NULL;
	}

	// larger buffer noticeably helps when reading gzip compressed genomes
	gzbuffer(file->fp, 1u << 20);

	file->kseq = kseq_init(file->fp);
	if (file->kseq == NULL) {
		log_error("cannot initialize the parser for %s", path);
		gzclose(file->fp);
		free(file);
		return NULL;
	}

	file->path = camil_strdup(path);
	return file;
}

int seq_read(struct seqfile *file, struct seqrec *rec) {
	int ret = kseq_read(file->kseq);

	if (ret >= 0) {
		rec->name = file->kseq->name.s != NULL ? file->kseq->name.s : "";
		rec->seq = file->kseq->seq.s != NULL ? file->kseq->seq.s : "";
		rec->len = file->kseq->seq.l;
		return 1;
	}

	// -1 is a clean end of file, anything else signals malformed input
	if (ret == -1) {
		return 0;
	}

	log_error("malformed or truncated record in %s (kseq code %d)",
	          file->path != NULL ? file->path : "input", ret);
	return -1;
}

void seq_close(struct seqfile *file) {
	if (file == NULL) {
		return;
	}
	if (file->kseq != NULL) {
		kseq_destroy(file->kseq);
	}
	if (file->fp != NULL) {
		gzclose(file->fp);
	}
	free(file->path);
	free(file);
}

int seq_readable(const char *path) {
	FILE *fp = fopen(path, "rb");

	if (fp == NULL) {
		return 0;
	}
	fclose(fp);
	return 1;
}

char *seq_basename(const char *path) {
	const char *begin;
	const char *slash;
	char *name;
	size_t len;
	size_t i;

	slash = strrchr(path, '/');
	begin = slash != NULL ? slash + 1 : path;

	len = strlen(begin);
	name = (char *)malloc(len + 1);
	if (name == NULL) {
		return NULL;
	}
	memcpy(name, begin, len + 1);

	// strip a compression suffix first so that "x.fa.gz" loses both parts
	if (len > 3 && strcmp(name + len - 3, ".gz") == 0) {
		len -= 3;
		name[len] = '\0';
	}

	// then strip the final extension, if any is left
	for (i = len; i > 0; i--) {
		if (name[i - 1] == '.') {
			name[i - 1] = '\0';
			break;
		}
	}

	// never return an empty name; fall back to the untouched basename
	if (name[0] == '\0') {
		free(name);
		return camil_strdup(begin);
	}

	return name;
}
