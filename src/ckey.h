#ifndef CAMIL_CKEY_H
#define CAMIL_CKEY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "camil.h"
#include <stdint.h>


// Derivation of the 64 bit key a core is stored and looked up under.
//
//   high 32 bits   core.label, the value lcptools already computed
//   low  32 bits   a hash of the core's own bases
//
// The label alone is a 32 bit hash, and an index of a few hundred million
// cores fills enough of that space for accidental matches to distort a run.
// Hashing the bases the core actually spans widens the key to 64 bits, where
// two distinct cores collide with probability around 1e-11. Put differently:
// with this key, equal keys mean equal substrings, up to that probability.
//
// Two details make the difference between this working and silently not:
//
// Reverse complement coordinates. init_lps2() produces cores identical to
// parsing the reverse complement of the sequence, and their start/end are
// positions in *that* string, not in the one that was passed in. Position j
// there is seqlen-1-j here, complemented. So an RC core hashes the complement
// of seq[seqlen-end .. seqlen-start) walked backwards, which is exactly the
// stream the parser consumed. Hashing seq[start..end) instead would produce
// keys that never match anything, with no other symptom.
//
// Case. LCP_INIT maps 'a' and 'A' alike, so the parser does not care, but a
// hash over raw bytes would: a soft masked reference such as GRCh38 stores
// repeats in lower case, and reads arrive in upper case. The hash therefore
// runs over the two bit alphabet codes rather than over ASCII, which also
// gives anything outside ACGT one fixed code instead of a negative one.
//
// The margin extends the hashed span by that many bases on each side, so a
// core only matches when its flanks match too. It raises specificity at the
// cost of sensitivity, since a difference anywhere in the flank breaks the
// match and cores closer than the margin to the end of a read or a contig
// cannot be keyed at all. It also makes repeated elements stop collapsing
// into one entry, so the index grows. Default 0.

// code for anything that is not A, C, G or T
#define CAMIL_CODE_OTHER 4

// Fills the lookup tables. Call once, next to LCP_INIT(), before any key is
// computed.
void camil_key_init(void);

// forward and complementing base codes, indexed by raw character
extern uint8_t camil_code[256];
extern uint8_t camil_code_rc[256];

// Computes the key of one core.
//
// `seq` and `seqlen` are the sequence the core was parsed from, `rc` says
// whether it came from init_lps2(), and `margin` is the number of flanking
// bases to fold in. Returns 0 on success, or -1 when the margin does not fit
// inside the sequence, in which case the core has no comparable key and the
// caller should skip it.
int camil_core_key(const char *seq, uint32_t seqlen, const struct core *cr, int rc, uint32_t margin, uint64_t *key);

#ifdef __cplusplus
}
#endif

#endif
