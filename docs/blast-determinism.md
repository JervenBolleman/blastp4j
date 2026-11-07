# How deterministic is BLAST+ blastp?

blastp4j aims for output identical to BLAST+ 2.17.0. That raises the question of how stable BLAST+'s own output is: across
platforms, thread counts and query batches. This note records what was measured on the bean data set (the first 500
proteins of the *Phaseolus vulgaris* proteome against the whole proteome, 28,388 sequences; `-outfmt '6 std btop'`,
55,629 lines), in November 2026.

In short: BLAST+'s output is deterministic for a given binary, platform, query file and batch size. It does not depend
on floating-point details. It does depend on three things that look harmless:

| Change | Effect on BLAST+'s output | blastp4j |
|---|---|---|
| Platform: macOS instead of Linux | 2 of 500 queries, from `qsort`'s order of ties | Linux's output on every platform |
| `-num_threads` | none on the 500 queries; 18 lines when one query fills its heap alone | does not depend on it |
| Query batch (`BATCH_SIZE`, which queries share a batch) | up to 10 of 500 queries; gained hits | follows the default batching only |

## Platform: `qsort` and ties

BLAST+'s macOS build (aarch64-macosx) differs from its Linux builds (x86-64 and aarch64, which agree with each other)
on 2 queries, found by the macOS workflow (`mac-compare.yml`, run 42301834031, Apple M1):

* `V7AVB7` against `V7D0Z4`: a different alignment, 59-437 / 349-734 (112 bits) on Linux, 83-437 / 347-734 (120 bits) on
  macOS.
* `V7BQ25` against `V7CQY7`: a different alignment (211 bits on Linux, 202 on macOS), which changes the order of 7 of
  that query's hits. Against `V7BSJ2`: the same coordinates and score, a different path (BTOP).

blastp4j on the same M1 gives Linux BLAST+'s output byte for byte, so the difference is in BLAST+.

**It is not floating point.** The composition-based stage rounds an adjusted score matrix for each query/subject pair,
the obvious place for a platform's floating point to show. Of 100,366 such matrices in this search, the score closest
to a rounding boundary (x.5) is 8.2e-9 away; 2 are within 1e-8, 75 within 1e-6. The two pairs that differ on macOS are
0.0027 and 0.0003 away, far beyond any difference in the last bits of `exp` or `log`.

**It is `qsort`.** glibc's `qsort` is a merge sort and keeps elements that compare equal in their input order. Apple's
is a quicksort (FreeBSD's, Bentley and McIlroy), which does not. `tools/qsort-order/apple-qsort.c` reimplements Apple's
algorithm; preloaded into the Linux x86-64 BLAST+ binary it reproduces the macOS output byte for byte, all 55,629 lines:

    gcc -O2 -shared -fPIC -o apple-qsort.so tools/qsort-order/apple-qsort.c
    LD_PRELOAD=$PWD/apple-qsort.so blastp -query queries.fasta -db proteome -outfmt '6 std btop' -num_threads 3

Running both sorts on every call (glibc's result kept) shows which call sites see ties in a different order: of 11
`qsort` call sites in this search, one. It is the sort by start in `Blast_HSPListPurgeHSPsWithCommonEndpoints`, which
orders HSPs by query and subject start, score and ends. Of its 182,212 calls, 497 come out in a different order. HSPs
equal on all those keys can still differ in their path or in the seed that the traceback stage starts from, and the
purge keeps only the first.

**How much depends on it.** `tools/qsort-order/reverse-ties.c` keeps glibc's order but reverses every run of equal
elements, the opposite choice for every tie. With it, 36 of the 500 queries change (80 query/subject pairs, 161
lines): 62 alignments keep their coordinates and score but take another path, 14 change score (by -34 to +8 bits), 3
keep their score but move, and one subject gets a different number of HSPs.

blastp4j sorts with `List.sort`, which is stable like glibc's merge sort, so it gives Linux's order on any platform.
Other C libraries are likely to give other results too: musl's `qsort` (Alpine) is a smoothsort and Windows' is not
stable either. Neither was tested.

## Threads

On Linux x86-64, `-num_threads` 1, 2, 4 and 8 (8 twice) and `-mt_mode 1` with 8 threads all give the same output for
the 500 queries, byte for byte.

Threads do matter when a query fills its 500-subject heap and is alone in its batch. Searched alone, `V7C040` (1,150
significant subjects) gives the same output with 1 thread and with `-mt_mode 1` (threads split the queries), but 18
lines differ with 2, 4 or 8 threads in the default `-mt_mode 0` (threads split the database). Each thread then keeps its
own heaps, and the composition-based stage skips different weak subjects once a heap is full. blastp4j follows the
single-threaded result whatever `-num_threads` is.

## Query batches

blastp concatenates its queries into batches of at least 10,000 residues (`BATCH_SIZE` in the environment overrides
this) and searches each batch as one sequence. All with 1 thread:

| Batching | Queries that differ from the default (of 500) |
|---|---:|
| `BATCH_SIZE=1000` or `4000` | 0 |
| queries in reverse order | 0 (the lines come in another order) |
| one query per run | 1 (`V7C040`, the heap effect above) |
| `BATCH_SIZE=30000` | 4 |
| `BATCH_SIZE=1000000` (one batch) | 10 (310 lines; most gain a hit) |

The batch effect is not in seeding: with `-ungapped` the hits of the affected query are the same in every batch. It is in
the gapped stage, and it depends on the exact contents of the batch rather than its size. Without composition-based
statistics (`-comp_based_stats 0 -evalue 1000`), `V7BHX8` gets 457 hits alone and in batches of queries 1-60, 1-72,
1-74, 1-100, 2-73 or 56-73, but 495 in the batch of queries 1-73. Adding query 73 (`V7C6L7`, 900 residues) to queries
1-72 changes three other queries: `V7B907` gains 89 lines, `V7BHX8` 38, and `V7B1R6`, whose 500-subject list is full,
exchanges some subjects (27 lines change). With
composition-based statistics this shows as one more hit for `V7BHX8` (`V7CPJ1`, E-value 7.3) with `BATCH_SIZE=30000`.

The likely cause, from reading the source but not verified by instrumenting BLAST+: the interval tree that decides which
seeds still need a gapped alignment (`blast_itree.c`) spans the whole concatenated batch, so its midpoints move whenever
a query is added or removed. Its containment test is exact, but the test for HSPs with a common endpoint, which removes
the worse of the two from the tree, is not: it follows one path down, stops at an HSP whose endpoint lies exactly on a
midpoint, and only looks at the first leaf it reaches. BLAST's own comment says the worst outcome is "a rare extra
gapped alignment", and the affected queries mostly gain hits. blastp4j's `HspTree` is a list per query that finds every
common endpoint, which matches BLAST+ with the default batches on all data sets so far. With `BATCH_SIZE=30000` blastp4j
still gives the default output, so it does not reproduce BLAST+'s result for that batch size.

## Consequences for comparisons

* The reference for blastp4j is Linux BLAST+ (glibc). Comparing with BLAST+ on macOS reports the 2 queries above as
  different, as the macOS workflow does.
* Single-threaded BLAST+ (or `-mt_mode 1`) is the reference for queries that fill their heap alone.
* Keep the default batch size. `BATCH_SIZE` is honoured only for the heap skipping that blastp4j models.
