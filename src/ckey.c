#include "ckey.h"


uint8_t camil_code[256];
uint8_t camil_code_rc[256];

// FNV-1a constants, used to fold the base codes into an accumulator one at a
// time. The stream is short, a few dozen codes, so a byte at a time loop is
// not worth vectorizing; the accumulator is finalized with a stronger mix.
#define CKEY_SEED  0xCBF29CE484222325ULL
#define CKEY_PRIME 0x100000001B3ULL

static inline uint64_t ckey_step(uint64_t hash, uint8_t code) {
	return (hash ^ (uint64_t)code) * CKEY_PRIME;
}

void camil_key_init(void) {
	int i;

	for (i = 0; i < 256; i++) {
		camil_code[i] = CAMIL_CODE_OTHER;
		camil_code_rc[i] = CAMIL_CODE_OTHER;
	}

	camil_code[(unsigned char)'A'] = camil_code[(unsigned char)'a'] = 0;
	camil_code[(unsigned char)'C'] = camil_code[(unsigned char)'c'] = 1;
	camil_code[(unsigned char)'G'] = camil_code[(unsigned char)'g'] = 2;
	camil_code[(unsigned char)'T'] = camil_code[(unsigned char)'t'] = 3;

	camil_code_rc[(unsigned char)'A'] = camil_code_rc[(unsigned char)'a'] = 3;
	camil_code_rc[(unsigned char)'C'] = camil_code_rc[(unsigned char)'c'] = 2;
	camil_code_rc[(unsigned char)'G'] = camil_code_rc[(unsigned char)'g'] = 1;
	camil_code_rc[(unsigned char)'T'] = camil_code_rc[(unsigned char)'t'] = 0;
}

int camil_core_key(const char *seq, uint32_t seqlen, const struct core *cr, int rc, uint32_t margin, uint64_t *key) {
	uint64_t hash = CKEY_SEED;
	uint32_t start = (uint32_t)cr->start;
	uint32_t end = (uint32_t)cr->end;
	uint32_t lo;
	uint32_t hi;
	uint32_t i;

	// The span is half open, and end-start is the number of bases the core
	// covers. Without room for the flanks the core cannot be keyed the same
	// way it would be elsewhere, so it is refused rather than keyed short.
	if (start < margin || end + margin > seqlen) {
		return -1;
	}
	lo = start - margin;
	hi = end + margin;

	if (!rc) {
		for (i = lo; i < hi; i++) {
			hash = ckey_step(hash, camil_code[(unsigned char)seq[i]]);
		}
	} else {
		// Walk the forward sequence backwards, complementing, which reproduces
		// the reverse complement stream the parser read.
		for (i = seqlen - lo; i > seqlen - hi; i--) {
			hash = ckey_step(hash, camil_code_rc[(unsigned char)seq[i - 1]]);
		}
	}

	*key = ((uint64_t)cr->label << 32) | (uint32_t)(camil_mix64(hash) >> 32);
	return 0;
}
