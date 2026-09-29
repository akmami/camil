# camil

Read classification with locally consistent parsing cores.

`camil` parses reference genomes into LCP cores with [lcptools](https://github.com/BilkentCompGen/lcptools), records every core of every reference once in a full index, carves out of it the cores that discriminate between whichever species a study is about, and assigns sequencing reads to the genome that owns most of their cores. 
It is the generalization of a three genome prototype (`misc/lammiq.cpp`) into a tool that takes any number of references, saves and reuses its index, and scales across threads.

## How it works

**Indexing.** Every reference is read sequence by sequence and parsed into LCP cores at a fixed level. 
Each core is reduced to its 64 bit lcptools label and recorded against the species it came from; a core seen a thousand times in one genome is stored once. 
Unless `--no-rc` is given, each reference sequence is also parsed as its reverse complement so that reads from either strand match. 
Nothing is filtered at this stage: the *full index* holds every core with the complete list of species it occurs in, so a collection of references is parsed once and reused.

**Subsetting.** A classify index is extracted from the full index for the species a run is about. 
A core is kept when it occurs in at least one of the chosen species and in at most `--max-share` of them; species that were not chosen are ignored entirely, so a core that is unique among the chosen species counts as unique even if the full index shares it with a hundred others. 
At the default of 1 only cores unique to a single chosen reference survive, which is the strict rule the prototype used. 
Raising it keeps cores shared by a few references, trading specificity for sensitivity. 
There is no ceiling: a shared core is a run of consecutive entries rather than a fixed size list.

**Classification.** A read is parsed at the same LCP level and every core is looked up. 
Each matched core votes for the species it names. 
The read goes to the genome with the strictly highest vote count, provided it clears `--min-hits` and `--min-ratio`; ties are never broken. 
Reads with no matched core at all are reported separately from reads whose cores disagree, so "nothing in the index looks like this" stays distinct from "several references look like this".

## Building

Requires a C99 compiler, zlib and pthreads. lcptools is a submodule and is built into a private prefix, so nothing is installed system wide.

```sh
git clone --recursive https://github.com/akmami/camil.git
cd camil
make          # builds deps/lcptools, then camil
make test     # end to end test on synthetic data
```

If the repository was cloned without `--recursive`, `make` initializes the submodule itself. 
lcptools is tracked on its `biolcp` branch and built as the variant `LABEL=64 POS=32 DCT=1 CORE=var ALPHABET=dna`; after pulling a new lcptools commit, run `make deep-clean && make` so that the library and its `config.h` are regenerated together.

Useful knobs:

```sh
make -j8                 # parallel build
make debug               # -O0 -g with the address and UB sanitizers
make install PREFIX=~/.local
make CFLAGS="-O2 -g"        # override the compiler flags
```

## Usage

```
camil index     -o <full index> [options] <genome.fa[,name]> ...
camil subset    -i <full index> -o <index> [options] <name> ...
camil classify  -i <index> [options] <reads.fq> [reads.fq ...]
camil run       -g <genome.fa[,name]> [-g ...] [options] <reads.fq> ...
```

Build the full index once, then extract a classify index for the species of interest and reuse it:

```sh
camil index -o refs.cfull -l 4 -t 16 alpha.fa beta.fa gamma.fa delta.fa
camil subset -i refs.cfull -o abg.cidx -t 16 alpha beta gamma
camil classify -i abg.cidx -o reads.tsv -s summary.tsv -t 16 reads.fq.gz
```

Species are named in `subset` by the short names the full index reports them under (see below), inline or one per line in a file given with `-S`. 
The output numbers them 0, 1, ... in the order they were named, which is also the row order of the summary.

Or build a classify index in memory and classify in one process, which is what the prototype did:

```sh
camil run -l 4 -t 16 -g alpha.fa -g beta.fa -g gamma.fa reads.fq.gz
```

Inputs may be FASTA or FASTQ, plain or gzip compressed; the format is detected from the file.

### Genome names

Reports label each genome with a short name. By default that is the file name with the directory and extension removed, which is unreadable once the paths are long, so a name can be given after a comma:

```sh
camil run -g /data/refs/GRCh38.p14.fa.gz,human -g /data/refs/GRCm39.fa.gz,mouse reads.fq.gz
camil index -o refs.cidx /data/refs/GRCh38.p14.fa.gz,human /data/refs/GRCm39.fa.gz,mouse
```

With many references, list them in a tab separated file instead, one `path` per line optionally followed by a tab and a name. Blank lines and lines starting with `#` are ignored:

```
# refs.tsv
/data/refs/GRCh38.p14.fa.gz	human
/data/refs/GRCm39.fa.gz	mouse
/data/refs/GRCz11.fa.gz
```

```sh
camil run -G refs.tsv reads.fq.gz
```

`-g`, `-G` and positional genomes can be mixed; genome ids follow the order they were given in. camil logs the mapping as it indexes, so the report can always be traced back to a file:

```
[   0.000] [INFO] genome 0: human (/data/refs/GRCh38.p14.fa.gz)
[   0.000] [INFO] genome 1: mouse (/data/refs/GRCm39.fa.gz)
[ 412.331] [INFO] indexed human (/data/refs/GRCh38.p14.fa.gz): 24 sequences, 3117275501 bases, 219402118 cores
```

Names must stay unique or two rows of a summary would be indistinguishable. The first genome to claim a name keeps it; a later one falls back to its file name, and if that is taken too, to a numbered variant. Every fallback is warned about. The split for `-g` is at the last comma, so a path that contains a comma has to go through `-G`, where the separator is a tab.

### Options

| Option | Applies to | Meaning |
| --- | --- | --- |
| `-g, --genome FILE[,NAME]` | run | reference genome, repeat once per genome |
| `-G, --genome-list FILE` | index, run | tab separated `path [name]` list of genomes |
| `-S, --species-list FILE` | subset | species names, one per line, first tab separated field |
| `-l, --level INT` | index, run | LCP level for references and reads (default 4) |
| `-n, --max-share INT` | subset, run | chosen species a core may occur in (default 1) |
| `-t, --threads INT` | all | worker threads (default 1) |
| `--no-rc` | index, run | do not index reverse complements |
| `-o, --output FILE` | index | where to write the full index |
| `-o, --output FILE` | subset | where to write the classify index |
| `-o, --output FILE` | classify, run | per read report |
| `-s, --summary FILE` | classify, run | summary table (default stdout) |
| `-i, --index FILE` | subset | full index to extract from |
| `-i, --index FILE` | classify | classify index to load |
| `--save-index FILE` | run | also write the index |
| `--min-hits INT` | classify, run | cores the winner needs (default 1) |
| `--min-ratio FLOAT` | classify, run | share of matched cores the winner must hold (default 1.0) |
| `-v, --verbose` | all | print progress and debug messages; without it only warnings and errors are printed |

### Output

The summary is a tab separated table, one block per read file:

```
# reads.fq.gz
#category	reads	percent
human	412033	41.20
mouse	389114	38.91
zebrafish	151902	15.19
ambiguous	31877	3.19
unclassified	15074	1.51
total	1000000	100.00
```

The per read report has one row per read, in input order:

```
#read	length	cores	matched	best_hits	second_hits	assignment	status
read_00001	151	6	5	5	0	human	assigned
read_00002	151	6	0	0	0	-	unclassified
read_00003	151	7	4	2	2	-	ambiguous
```

## Design notes

### The key

Every core is stored under the 64 bit label lcptools computes for it: at level 1 the packed bases of the core, and at every level above a 64 bit MurmurHash3 over the labels of the cores it was compressed from and their count. 
lcptools has to be built with `LABEL=64` for that; the Makefile pins the variant (`LCP_VARIANT`) and `camil.h` refuses to compile against anything else, so a 32 bit label cannot silently be zero extended into a 64 bit key.

`LCP_INIT` maps `a` and `A` alike, so a soft masked reference such as GRCh38 matches uppercase reads. 
`init_lps2` produces cores identical to parsing the reverse complement, and its labels are compared on equal footing with the forward ones, so reads from either strand match.

### Building

Collection is lock free. A worker parsing a sequence appends one twelve byte record per core to its own sink, which keeps a separate chunk list per bucket, where the bucket is the top bits of the key. No shared counter, no atomic, no shared cache line.

Reduction is lock free too. Each bucket is owned by one thread, which gathers that bucket's chunks from every sink, radix sorts them and collapses repeats of a core within one species. Buckets are disjoint key ranges, so the sorted buckets concatenated in order are globally sorted: no merge, no heap. For a classify index built by `run`, cores present in more than `--max-share` species are dropped here; for the full index nothing is.

### The full index

The full index is a hash table from core key to a run in one flat array of species ids: `key -> (offset, count)`, the species of the core being `sids[offset .. offset+count)` in ascending order. 
The table is a [khashl](https://github.com/attractivechaos/klib) ensemble, split into sub tables by the low bits of the key hash so that each stays inside khashl's 32 bit addressing while the whole may hold billions of cores. 
The value packs the offset into 40 bits and the count into 24, which caps an index at a trillion (core, species) pairs and 16 million species. 
On disk the sub tables are dumped as they sit in memory, bitmap and buckets, so loading is a read with no rehashing.

Extracting a subset is one pass over the table: for each core, count how many of its species were chosen, keep it when that is between 1 and `--max-share`, and push the chosen (key, species) pairs through the same parallel radix sort a fresh build uses. 
That turns the hash order of the table into the sorted classify table without a separate sort step, and renumbers the species 0, 1, ... in the order they were named.

The cost of keeping everything is size: about 16 bytes of table per distinct core at a 75% load, plus 4 bytes per (core, species) pair, all of it held in memory while a subset is extracted. 
The worked example below, one human genome next to two salmonella strains, comes to roughly 5 GB in memory and on disk.

### Querying

The index is a sorted array of keys, a parallel array of values, and a directory that turns the top bits of a key into the short run of keys it could be in. Nothing is written after building, so lookups take no locks and the arrays stay shared in every core's cache.

Keys and values are separate arrays because most lookups miss, and a miss reads one cache line of keys and never touches values. A lookup is two dependent misses, the directory then the keys, so classification prefetches both levels several cores ahead — a read is reduced to all of its keys before any is looked up, which is what makes that possible. Measured against the previous open addressed table, the two stage prefetch hides the extra miss completely.

Classification is arranged so a worker writes only memory nobody else touches: reads are batched, each batch is cut into contiguous chunks aligned to whole cache lines of result records, and each worker keeps its own vote counters. Only the species a read actually hit are examined when scoring it, so the scan costs nothing when an index holds thousands of them. Two batches are kept in flight, so decompressing the next overlaps with classifying the current. Output is written in input order, so results do not depend on the thread count.

### Index size

The three arrays of a classify index are exactly as long as there is data: twelve bytes per core, plus about half a byte per core for the directory. Nothing is rounded to a power of two and there are no empty slots.

What drives the size is how many cores survive: roughly one distinct core per 14 bases of reference at level 4, counting both strands. The sharing limit only removes cores shared by more than `--max-share` chosen species, so it does nothing at all for references that are not related to each other. A human genome next to a few bacterial ones keeps essentially every human core, and the human contributes over 99% of the table.

Worked example, one human genome and two salmonella strains at level 4: ~220M cores, of which ~0.7M are salmonella, for a classify index around 2.8 GB.

To make one smaller: choose fewer species in the subset, raise `--level` (fewer, longer cores, but fewer cores per read), or `--no-rc` and accept that reads off the reverse strand will not match.

Peak memory while building is dominated by lcptools rather than by camil: `init_lps` allocates a core array of `len/1.5` records for the whole sequence, about 16 bytes per base, so a 250 Mbp chromosome costs roughly 4 GB while it is being parsed. `CAMIL_INDEX_INFLIGHT_BASES` caps how much sequence is queued at once and is what to lower if a machine runs short.

### Collisions

Unrelated cores can still collide, at `entries / 2^64`. For an index of 220M cores that is one in 1e11 per core probed, against roughly one in 240 for a 32 bit label. Measured with the earlier 64 bit key, 163,004 cores from unrelated reads produced 0 spurious matches where a 32 bit label predicts ~678.

Consequences worth knowing:

- at most 16,777,214 species per index, the width of the count field in the full index
- `--max-share` is unbounded, since a shared core is a run of entries rather than a fixed size list

## Layout

```
include/      public headers, one per subsystem, each documented at the top
src/          implementation
  cbuild.c      lock free collection, then a parallel radix sort per bucket
  cfull.c       the full index: khashl table plus flat species array, subset extraction
  ctable.c      sorted keys plus a directory, queried while classifying
  index.c       collection from genomes, classify index serialization
  classify.c    batching, threading, decision rule
  seqio.c       FASTA/FASTQ reader (the only user of kseq.h)
  tpool.c       fixed size thread pool
  opt.c         command line parsing
  logger.c      diagnostics and timing
deps/lcptools   submodule, branch biolcp
deps/klib       kseq.h, ksort.h and khashl.h, vendored from klib (MIT)
misc/lammiq.cpp   the original three genome prototype
tests/            end to end test
```

## Citing

If you use camil, please cite the LCP work it is built on:

> LCPan: efficient variation graph construction using locally consistent
> parsing. Akmuhammet Ashyralyyev, Zülal Bingöl, Begüm Filiz Öz, Salem Malikic,
> Uzi Vishkin, S. Cenk Sahinalp, Can Alkan. Genome Biol (2026).
> https://doi.org/10.1186/s13059-026-04088-w
