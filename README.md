# WARNING

Warning on this *SLOP* this is slop driven development. Do *NOT* use it for production or anything.
This is turning into an investigation on register use on x86 by C2 and Graal compared to GCC that has little to do with correct blast results (other than as an experimental constant).

# blastp4j
A protein database search in Java that aims to reproduce the output of NCBI BLAST+ `blastp`: the same hits, raw
scores, bit scores, E-values and alignments, for the same parameters. It builds on
[cssw4j](../cssw4j) (a Java Vector API port of the SSW Smith-Waterman library) and takes ideas from
[BLVector](https://pmc.ncbi.nlm.nih.gov/articles/PMC7884812/).

Status: the blastp search is implemented and reproduces blastp's tabular output exactly, with blastp's default
composition-based statistics (mode 2) and without (`-comp_based_stats 0`). For the 24 reference queries against
Swiss-Prot, all lines are identical (11,263 and 15,107), on every column including the alignment (BTOP), and in the
same order. The `blastp4j` command takes blastp's options and writes its tabular output.

## Why this is more than "SSW plus statistics"

To report the same alignments and scores as blastp, blastp4j has to run the same algorithm, not an exact
Smith-Waterman search: BLAST's gapped extensions are X-drop extensions from seeds, so their scores and end points can
differ from the optimal local alignment, and its default scores are adjusted per subject by composition-based
statistics. The pipeline to reproduce (task `blastp`, defaults):

1. **Lookup table.** All words of 3 residues whose BLOSUM62 score against a query word is at least 11 (the
   neighbourhood words), mapped to query positions.
2. **Scan and two-hit trigger.** Each subject position looks up its word; an extension starts when two non-overlapping
   hits lie on the same diagonal within 40 residues.
3. **Ungapped X-drop extension** (7 bits). HSPs reaching the gap trigger score continue.
4. **Preliminary gapped extension** (X-drop 15 bits, score only) from a start point inside the HSP.
5. **Composition-based statistics** (mode 2, Yu and Altschul 2005): per subject, decide whether to adjust the matrix
   to the query and subject composition, and redo the alignments with the adjusted matrix.
6. **Final gapped extension with traceback** (X-drop 25 bits).
7. **Statistics and hit list.** Bit scores from the Karlin-Altschul parameters of gapped BLOSUM62 11/1; E-values with
   Spouge's finite-size correction (per subject, from both sequence lengths and the Gumbel parameters, scaled to the
   database), which BLAST+ uses for BLOSUM62 instead of the classic length-adjusted `K m' n' e^(-lambda S)`; E-value
   cutoff 10, at most 500 subjects, BLAST's ordering.

The algorithms are described in the BLAST papers (Altschul et al. 1990, 1997; Altschul et al. 2005; Yu and Altschul
2005) and implemented in the BLAST+ source code, which is a US government work in the public domain, so it can be
followed closely.

## Where cssw4j and BLVector fit

* cssw4j's striped Smith-Waterman gives an **exhaustive baseline**: the best local alignment of the query against every
  subject. It shows how many real hits the heuristic misses and is the reference for "how good can it get".
* BLVector replaces BLAST's lookup table and two-hit trigger by a brute-force, vectorized 4-mer filter followed by full
  Smith-Waterman. That gives different hits from BLAST, so it does not fit "match BLAST+ closely" as the main path, but
  its ideas (vectorized seed scanning, a reduced alphabet) can be tried as optional, measured alternatives. BLVector is
  AGPL-3.0 licensed C code; only its published ideas are used here, none of its code.

## Plan

Each milestone is checked against the BLAST+ reference results, and each experiment gets its own JMH benchmark and
commit.

1. **Statistics** (done). Bit score and Spouge E-value for BLOSUM62 11/1, and BLAST's tabular number formatting.
   Check: the `bitscore` and `evalue` columns of all 15,107 hits of the `-comp_based_stats 0` reference are recomputed
   exactly from the raw `score`, `qlen` and `slen` columns. Cost: 0.3 us per alignment (0.7 us with formatting).
2. **Database and queries** (done). FASTA (plain or gzip) into memory in the layout of BLAST's `.psq` file: NCBIstdaa
   codes with a sentinel between sequences; ids are the first word of the header, as BLAST prints them. Check:
   Swiss-Prot reads as 575,748 sequences and 209,017,843 residues, byte-identical to `sprot.psq`; every subject id
   and length of the reference matches, and so do the 24 query ids and lengths. Loading Swiss-Prot takes about
   1.1 s from plain FASTA and 2.3 s from gzip.
3. **Exhaustive baseline** (done) with cssw4j over all of Swiss-Prot; see [the results](#exhaustive-baseline).
4. **Seeding** (done): lookup table, two-hit trigger, ungapped X-drop extension; see [the results](#seeding).
5. **Gapped X-drop extension and traceback** (done); see [the results](#gapped-extension-and-traceback).
6. **Composition-based statistics** (done, mode 2); see [the results](#composition-based-statistics).
7. **Command line** (done): `blastp4j`, with blastp's options and `-outfmt 6`; see [usage](#usage).
8. **BLAST databases** (done): `-db` also takes a database made by `makeblastdb`; see
   [BLAST databases](#blast-databases).

### Long-term to do: how far can almostblastp be from blastp?

almostblastp's changes (`util.Arithmetic.FAST`) change floating-point results in their last bits. So far its output
has been blastp's on every workload measured, but no bound says how often or by how much it can differ. A numerical
analysis should bound the combined difference of all FAST changes, from the floating-point error up to the output:

1. **The changes' own errors.** For each FAST change, the error per operation against BLAST's: the Vector API's
   logarithm (SVML on x86-64, SLEEF on AArch64; their documented ulp bounds, to check against `Math.log` over the
   inputs that occur), sums in another order (the bound for a reordered sum of n terms), fused multiply-adds (one
   rounding instead of two).
2. **Through the target-frequency optimizer.** The optimizer is Newton's method on a smooth convex problem, stopped
   when the residual is below 1e-8 (`COMPO_ADJUST_ERR_TOLERANCE`). A perturbation of its inputs and of each
   iteration moves the solution by about the condition number of the Newton system times the perturbation; the
   stopping rule itself allows a difference of about the tolerance between any two runs, which is likely larger than
   the rounding. Bound the difference in the target frequencies, then in the frequency ratios, then in the scores
   before rounding (log(ratio) / λ: the relative error of a ratio becomes an absolute error of the score, scaled by
   32 / λ). A FAST run can also take one Newton iteration more or fewer, or fail to converge where BLAST's does (then
   the matrix is not adjusted at all): measure how close to the tolerance the last residuals are.
3. **To the integer matrix.** Each adjusted score is rounded to an integer (scaled by 32). A difference ε before
   rounding flips the integer only when the value is within ε of a half-integer. Estimate the distribution of the
   distance to the nearest half-integer over real compositions (it should be about uniform), which gives the expected
   number of flipped entries per adjusted matrix, about 2ε for each of the 400 entries.
4. **To the alignments.** A flipped entry changes the score of every aligned pair of those two letters by ±1 (a 32nd
   of a raw score unit): the alignment score changes by at most the number of such pairs, usually 0 or 1. Then the
   E-value (from the scaled score), the order of hits with nearly the same score, the hits at the E-value cutoff and
   at `-max_target_seqs`, and the X-drop end points (a changed cell can let an extension stop or go on) can change.
   Bound each, and how often each happens given step 3.
5. **Check it empirically.** Run almostblastp against blastp4j on large workloads (all of Swiss-Prot against itself,
   the bean proteome, the 24 human queries): count differing adjusted matrices, flipped entries, and differing output
   lines, and compare them with the bounds. A test could then accept almostblastp's output when it is within those
   bounds of blastp's, instead of identical (`AlmostBlastpTest`).

## Usage

    mvn package -DskipTests
    bin/blastp4j -query queries.fasta -db blastdb/sprot -outfmt "6 std btop" -num_threads 8 -out hits.tsv
    bin/blastp4j -query queries.fasta -db uniprot_sprot.fasta -outfmt "6 std btop" -num_threads 8 -out hits.tsv

* The database is a BLAST protein database (the name given to `makeblastdb -out`), or a FASTA file, plain or
  gzip-compressed. A BLAST database is memory-mapped, so it is not limited in size by the Java heap or arrays, and the
  operating system pages it in. A FASTA file is read onto the heap: Swiss-Prot (209 million residues) takes 2.3 s and
  0.9 GB, at most 2 GB of residues.
* From a FASTA file the results are those of blastp against a database made from it with `makeblastdb -parse_seqids`.
* Options, with blastp's names and defaults:
  * `-query`, `-db`, `-out`
  * `-outfmt "6 [fields]"`: the fields are `std qseqid sseqid qacc qaccver sacc saccver qlen slen qstart qend sstart
    send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovhsp btop
    stitle`
  * `-evalue` (10), `-max_target_seqs` (500), `-max_hsps`
  * `-comp_based_stats` (2; or 0)
  * `-num_threads`: the default is the number of processors. The output does not depend on it.
* `-task blastp`, `-matrix BLOSUM62`, `-gapopen 11`, `-gapextend 1`, `-word_size 3`, `-threshold 11`,
  `-window_size 40` and `-seg no` are accepted with these values only. Anything else is rejected with a message.
* Use a HotSpot JDK. `bin/blastp4j` uses `$JAVA_HOME` and passes this machine's options from `bin/jvm.options`, then
  `$BLASTP4J_OPTS`, to the JVM.
* `bin/almostblastp` (`cli.AlmostBlastp`) takes the same options. It is for changes that make the search faster but may
  change its results slightly, such as fused multiply-adds in the floating-point code; `blastp4j` stays exact. For now
  it is a copy of `cli.Blastp` (without the daemon) and gives the same results, which `AlmostBlastpTest` checks.
  The architecture comparison workflow times it in the full search next to blastp4j and BLAST+ (its `almostblastp`
  input, default true), and the report plots it against BLAST+.
  * Its first change (121ec20): the composition optimizer's logarithms (400 per Newton iteration) a vector at a time,
    with the Vector API's `LOG` (SVML on x86-64, SLEEF on AArch64), within about an ulp of `Math.log`, and the sums
    in another order (`util.Arithmetic.FAST`; blastp4j keeps `BLAST`). The results are still blastp's on the bean and
    human queries. On the runners (run 42096646036, against blastp4j in the same rounds): 1.2% less wall time on
    AArch64 and 1.4% on a Xeon 8573C, two of three rounds each; locally 1.3% fewer cycles on the bean queries.
  * *Tried, no gain:* the optimizer's Cholesky factorization and solve with fused multiply-adds and four partial
    sums (39f9f3e, reverted).
  * The optimizer a class of its own, `compo.OptimizeTargetFreqFast` (195c589): every loop of a Newton iteration a
    vector at a time (logarithms, gradients, the products with A and A^T, the residuals, the Cholesky factorization's
    and solve's dot products and updates, the step bound), and `Blast_FreqRatioToScore`'s logarithms too (3e10055).
    BLAST's `OptimizeTargetFreq` keeps no trace of it; `OptimizeTargetFreq.of` picks the class. The results are still
    blastp's on the bean and human queries. On the runners (run 42114578082): 2.7% less wall time than blastp4j on
    AArch64 and 4.0% on an EPYC 9V74, every round; locally 4.4% fewer cycles on the bean queries.
  * *Tried for completeness, no gain:* BLOSUM62's logarithms as constants in the fast optimizer (03af68a, on the
    branch `exp/logq`, not merged): a hidden copy whose class data is log q, so that the starting scores take 40
    logarithms instead of 400 (log q − log row − log column) and each iteration's log(x / q) is log x − log q. Run
    42126697919: almostblastp over blastp4j 0.969 on AArch64 and 0.957 on an EPYC 9V74, against 0.970 and 0.963 in the
    run before (another x86-64 CPU): within the noise, as expected (the setup is about a seventh of the optimizer's
    logarithms). Output identical.
  * `-XX:+TrustFinalNonStaticFields`, which lets C2 treat every final instance field as a constant (by default only
    those of records, hidden classes and some JDK classes), was tried (9ee5ac4, reverted): within the noise on the
    GitHub runners (1.4% less wall time on AArch64, 0.2% on x86-64) and here.
* `bin/fastblastp` (`cli.FastBlastp`) gives the results of BLAST+'s `blastp -task blastp-fast`, identical to it, with
  the same options; `-task blastp-fast`, `-word_size 5` and `-threshold 20` are accepted with these values only. A copy
  of `cli.AlmostBlastp` with exact arithmetic and `search.Task.BLASTP_FAST`, which changes two things:
  * The seeding: words of 5 letters of BLAST's 15-letter compressed alphabet (`seed.CompressedAlphabet`: the groups
    `ST IJV LM KR EQZ A G BD P N F Y H C W`, scored from BLOSUM62's frequency ratios), neighbours scoring at least 20
    (`seed.CompressedTable`), scanned by `seed.FastBatchWordFinder`; the two-hit logic is a copy for 5-letter words
    (`seed.FastWordFinder`). The table keeps BLAST's order of a word's hits (the first 3, then chunks of 4, the last
    first), which decides ties later.
  * Before the gapped extensions, BLAST's chaining estimate (`search.Chaining`, `s_ChainingAlignment`) drops the
    ungapped alignments that cannot reach the gapped cutoff.

  On the bean queries blastp-fast takes 4.4 times fewer cycles than blastp but misses 30% of its query-subject pairs
  (7.2% of those with E < 1e-10); see [docs/blast-determinism.md](docs/blast-determinism.md) for the comparison with
  4-letter words. fastblastp's output is BLAST+'s: the bean queries (30,538 lines; 39,448 with `-comp_based_stats 0`),
  the human queries against Swiss-Prot (8,706 lines), 2,500 bean queries, a single query, and the CLI sample
  (`FastBlastpTest`). The architecture comparison workflows time both in the full search (their `fastblastp` input,
  default true) and compare their outputs with each other; see [fastblastp's speed](#fastblastps-speed).
* The word finder uses the Vector API, still the incubator module `jdk.incubator.vector`: `bin/jvm.options` adds it,
  and the JVM then prints a warning about incubator modules on stderr. Other ways of running blastp4j need
  `--add-modules jdk.incubator.vector` too.

Checked against blastp:

| Run | Lines | Identical |
|---|---:|---:|
| 24 queries vs Swiss-Prot, defaults, the reference columns | 11,263 | all |
| 24 queries vs Swiss-Prot, `-comp_based_stats 0` | 15,107 | all |
| same two runs with `-db` the Swiss-Prot BLAST database instead of the FASTA file | 26,370 | all |
| 4 queries, all supported columns (`qseq`, `sseq`, `positive`, `ppos`, `qcovhsp`, ...) | 3,007 | all |
| 4 queries, default `-outfmt 6` (std: `qaccver saccver ...`) | 3,007 | all |
| 4 queries, `-evalue 1e-5 -max_target_seqs 50` | 397 | all |
| 4 queries, `-comp_based_stats 0 -evalue 100 -max_target_seqs 1000 -max_hsps 1` | 1,974 | all |
| 24 queries vs a 400-sequence database, `6 std btop` (`BlastpTest`) | 196 | all |
| same, `-comp_based_stats 0 -evalue 1000 -max_target_seqs 20`, other columns (`BlastpTest`) | 1,422 | all |
| same database in six volumes (`.pal`), `6 std btop` (`BlastpTest`) | 196 | all |
| 10 identifier kinds, `sseqid sacc saccver stitle`, with and without `-parse_seqids` (`BlastpTest`) | 20 | all |

Details that output depends on:

* `std` is `qaccver saccver ...`, not `qseqid sseqid ...`. With a `-parse_seqids` database the subject is printed as
  its accession (`P60052`). The query, read without `-parse_deflines`, keeps its whole identifier.
* `qcovhsp` is rounded below 99% and truncated above, so 684 of 685 residues is 99, not 100.
* The hit-list sizes follow `-max_target_seqs`: 2 times it plus 50 (at most 50 more without composition-based
  statistics, at least 1,050 with).
* An `-evalue` above 10 turns off the restricted alignment of the preliminary stage.

End-to-end time for the 24 queries against Swiss-Prot (i7-6700HQ, Temurin 25; blastp against its BLAST database,
blastp4j including reading the FASTA file):

| | blastp4j | blastp |
|---|---:|---:|
| 1 thread | 143 s | 80 s |
| 8 threads | 36 s | 19 s |

blastp4j is about 1.8 times slower over a whole run. The single-query benchmarks above (`BlastpSearchBench`) showed a
smaller gap, 10 to 30%. They repeat one warmed-up query, and the blastp times they are compared with include starting
the program, so they understated the difference.

On JDK 27 (OpenJDK 27+35, C2, through `bin/blastp4j`, three interleaved rounds), blastp4j took 1.80 times BLAST+'s
user-mode cycles on 8 threads (724.7 G against 402.7 G), the same as on Temurin 25.0.4 (727.2 G). JDK 27's compact
object headers, on by default, save about 1% (730.8 G without).

With the vector word finder and the gapped-stage changes since then (main at 7888777, three interleaved rounds),
blastp4j on the JVM took 1.53 times BLAST+'s cycles (615.9 G against 403.2 G). A native image built with
profile-guided optimisation (Oracle GraalVM 25.4.4, `bin/build-native-pgo.sh`, trained on 50 bean proteins) took 1.67
times (672.5 G, with 35% more instructions than the JVM; not profiled yet). In a native image the word finder reads
the presence table with scalar loads.

With the word finder's constants copy and the score-only extensions in vector lanes (main at 36cbdf8, JDK 27 C2, three
interleaved rounds, output identical), blastp4j took 1.43 times BLAST+'s cycles (570.7 to 579.2 G against 401.5 to
402.7 G) and 1.27 times its instructions (650.6 to 664.0 G against 516.7 G). Wall time was 24.9 to 25.6 s against 19.0
to 19.3 s, and peak memory 767 to 788 MB against 357 MB.

Most of what remained beyond the word finder was the JIT: all 24 queries form one BLAST batch, so the
composition-based stage started after the whole preliminary stage and ran mostly before C2 had compiled it (49 G
cycles cold, 11 G warm, BLAST+ 12 G), and C2 spent 20 s of compiler-thread time, 5.7 s of it on
`OptimizeTargetFreq.optimize` with all its helpers inlined. A thread that warms up the composition-based stage during
the preliminary one (a316748; once per JVM since 2279159, since once per batch it cost 14% more instructions on
500 bean queries, which form many batches) and the optimizer's helpers left out of line (a2a0da3, in `bin/jvm.options`) bring
blastp4j to 1.32 times BLAST+'s cycles (main at 97df58a, three interleaved rounds, output identical: 525.8 to 538.1 G
against 402.5 to 411.7 G) and 1.22 times its instructions (625.2 to 637.1 G). Wall time was 23.0 and 23.3 s against
19.0 and 20.1 s in the rounds without throttling. Where the first batch is short, as with the 500 bean queries of the
architecture comparison workflow, the warm-up costs 1.5 to 2.5% of cycles instead.

### fastblastp's speed

A blastp-fast search is about 4 times shorter than a blastp one, so on the 500 bean queries the JIT dominated it:
locally (i7-6700HQ, 4 threads, `-XX:+CITime`) C1 and C2 compiled for 15.9 s of the 59.5 s of user time (27%), and
the runners showed fastblastp at 1.2 to 1.8 times BLAST+'s wall time. With 2,500 queries compilation is 10% (20.3 s of
211.3 s), and fastblastp took 59.0 s against BLAST+'s 55.2 s (1.07 times), with the same user time (211.3 s against
212.6 s). So the workflows give blastp-fast its own queries, the first `fast_queries` (default 2,500) bean proteins,
and run the reference commit's fastblastp in the same job (`fastref`), as they do for blastp4j (0aecb1f). In that
run the two identical builds differed by about 1% (AArch64 44.72 s and 45.17 s, x86-64 58.40 s and 58.42 s).

A profile after the warm-up (async-profiler attached from second 25 to 50 of the 2,500-query search): the
composition-based stage 58% (the traceback lanes 27%, the target-frequency optimizer 25%), the compressed scan 18%,
the preliminary lanes 14%, building the compressed tables 2.3%, compilation 2%. Per thread, in wall time, the main
thread waited for the workers 82.7% of its time, and spent the rest building the compressed tables (7.5%), computing
the queries' ungapped Karlin-Altschul parameters (about 5%) and writing the output (1.7%). Meanwhile the workers
waited: BLAST+ kept 3.8 of 4 cores busy, fastblastp 3.25 (Neoverse-N2) to 3.35 (Xeon 8370C). Hence these experiments,
each measured on the runners against its parent (`fastref`), output identical to BLAST+ in every round:

| Commit | Change | AArch64 (Neoverse-N2) | x86-64 | |
|---|---|---:|---:|---|
| f8ba9e9 | a presence bit per cell of the compressed table (95 KB), tested before its 3 MB start array, as BLAST tests its pv array | | | kept; -0.6% cycles locally (measured on the runners before `fastref`) |
| 748e7d2 | blastp's 3-letter lookup tables built only for the composition stage's warm-up, not for every batch | -4.6% | -4.5% (EPYC 7763) | kept, every round |
| b4a6547 | the scan as a hidden copy with its window and mask as constants, as `BatchWordFinder`'s | +0.7%, rounds mixed | +1.8% (EPYC 9V45), rounds mixed | reverted (37521f5) |
| 2939f0c | the queries' Karlin-Altschul parameters computed once (the CLIs computed them for the skip check and `prelim` again), and the ideal lambda a constant | -2.2% | -2.4% (EPYC 7763) | kept, every round |
| abdf0b4 | the next query batch prepared (skip checks, cutoffs, lookup tables) on a thread of its own while this one is searched | -4.7% | -3.9% | kept; blastp4j -2.9% on AArch64 |

Changes are in median wall time. The look-ahead costs 3% more cycles on a 4-core laptop, where its thread competes
with the four workers; it pays in wall time. Since then (run 42421494383, main at 37521f5, median of three rounds):

| | AArch64 (Neoverse-N2) | x86-64 (EPYC 9V45) |
|---|---:|---:|
| fastblastp against BLAST+ blastp-fast (2,500 queries) | 39.12 s against 32.33 s (1.21) | 27.73 s against 38.64 s (0.72) |
| blastp4j against BLAST+ (500 queries) | 29.43 s against 28.66 s (1.03) | 26.06 s against 35.31 s (0.74) |

Before the experiments (0aecb1f) fastblastp was at 1.35 on the N2 and 0.96 on a Xeon 8370C. The x86-64 ratio depends
on the runner's CPU: 0.85 on an EPYC 7763 in 748e7d2's run.

On the N2 the profile of fastblastp (run 42462383540) then had the scan at 23.6% of the samples, 15.0% of them on the
two instructions after the loads of the cell's start and of its first entry: two dependent cache misses per looked-up
position, the start array being 3 MB and the N2's L2 1 MB per core. Two more experiments, the second on its own:

| Commit | Change | AArch64 (Neoverse-N2) | x86-64 | |
|---|---|---:|---:|---|
| 81dfb15 | blastp's AArch64 deferral of the finder calls (see [Seeding](#seeding)) in the fast scan | -0.1% | +0.2% (EPYC 7763, deferral off) | reverted (fd5ee9b) |
| 11cd4d9 | each cell's first entry next to its start (`CompressedTable.cells`, 6 MB instead of 3), so that a cell with one hit is one line, as BLAST's backbone cell holds its first hits | -8% (fast 37.86 s against 43.89 s for `fastref`, whose 7d151d5 cost about 6%) | -1.9% (EPYC 9V45) | kept; scan 19.2% of the samples |
| ad222c6 | the scan in two passes: first the positions whose cell has hits, then their hits, each touching `cells` of the position 8 hits ahead, a software prefetch (Ainsworth and Jones, CGO 2017) | -1.1%, one round of three lost; scan 13.3% of the samples | -0.8% (EPYC 7763) | kept |
| 1d211d2 | 16 hits ahead | 0%; scan 12.0% | -0.6% (Xeon 8573C) | kept |
| baf7301 | 4 hits ahead | -1.7%, round 3 lost; scan 16.8% | +3.3% (EPYC 9V45), rounds mixed | reverted |
| 683d08a | the compressed table's entries as chars when the query has fewer than 65,536 positions (half the bytes) | +1.9%, rerun -1.0%, rounds mixed | +0.7%, rerun 0% (EPYC 9V45) | reverted (9af8a64) |

Locally (i7-6700HQ, where 3 MB fit in L3) the inline entries were neutral in cycles, with 15% more LLC misses. After
them (run 42473840187, median of three rounds) fastblastp took 37.86 s against BLAST+ blastp-fast's 34.09 s on the
N2 (1.11), and 27.93 s against 37.54 s on an EPYC 9V45 (0.74); blastp4j 29.19 s against 26.28 s (1.11; that runner's
BLAST+ was faster than run 42421494383's). The prefetch took the scan from 19.2% to 12.0% of the samples on the N2 but
the wall time only about 1%: the out-of-order core already overlapped part of the misses, as Ainsworth and Jones
found for out-of-order cores (1.1 times, against 2 to 3 on in-order ones). Chars for the entries did not gain either.
The first run's ARM machine-code listing showed the hit loop at 0% and a tenth of the samples outside the printed
methods: with PrintAssembly its compiles were late and it ran in the interpreter. The rerun's listing was normal;
`bin/machcode.py` now prints the share of the samples that its listing covers.

### Daemon

`blastp4j daemon` is a search server: it opens its databases once, takes jobs over HTTP (the JDK's
`com.sun.net.httpserver`, no dependencies) and runs them one at a time, each on all threads.

    bin/blastp4j daemon -dir jobs -db sprot=blastdb/sprot -db bean=proteome.fasta -port 8080 -num_threads 8
    curl --data-binary @queries.fasta 'http://localhost:8080/jobs?db=sprot&outfmt=6+std+btop'   # prints the job id, 1
    curl -i http://localhost:8080/jobs/1                       # 202 while queued or running, then 200 and the output
    curl http://localhost:8080/databases                       # name, sequences, residues

* A job's parameters are blastp's options (`outfmt`, `evalue`, `max_target_seqs`, `max_hsps`, `comp_based_stats`
  and the fixed ones), checked as on the command line; a bad job is answered with 400.
* A job is stored as `input/<id>.fasta` and `input/<id>.options` in the job directory. When it is done the output is
  `output/<id>.tsv`, or the error `output/<id>.error`, and the input is deleted. Ids count up from one more than the
  largest id in either directory, and jobs left in `input` run again at startup.
* A job's queries are searched in blastp's batches, as `blastp4j -query` does, so the output is blastp's for the job's
  FASTA file. Jobs are not merged into larger batches: with composition-based statistics the result depends on the
  batch, so a job's output would then also depend on the other jobs.
* It listens on 127.0.0.1 unless given `-bind`, and has no authentication.
* The word finders specialized at run time (hidden classes with a query's window, X-drop, cutoff and diagonal mask as
  constants) are kept in an LRU cache of `-Dblastp4j.specializations` (256) entries, so a long-running daemon does not
  keep one for every cutoff it meets.

The 24 human queries against Swiss-Prot, posted three times to one daemon (8 threads): 22.5, 19.6 and 19.6 s, output
identical to blastp's each time. The first job is about as fast as the command line, which compiles as it goes; the
later ones run compiled code (BLAST+: 19.0 to 20.1 s).

## BLAST databases

`db.BlastDbReader` reads the protein databases `makeblastdb` makes, format version 4 or 5:

* `.pin`, the index: title, number of sequences, total length, and the offsets of each sequence and header. All
  integers are big-endian except the total length, which is a little-endian 64-bit number.
* `.psq`, the sequences: NCBIstdaa codes with a 0 byte between sequences. This is the layout blastp4j uses, so the
  file is used as is, with the offsets from `.pin` as sequence starts.
* `.phr`, the headers: one `Blast-def-line-set` per sequence in binary ASN.1 (BER, with indefinite lengths).
  `db.Asn1` is a small BER reader for it. A header is decoded only when a hit is printed.
* `.pal`, an alias file that lists the volumes of a large database (`DBLIST`). Alias files that restrict a database
  (`OIDLIST`, `GILIST`, `SEQIDLIST`, `TAXIDLIST`, ...) are rejected. The LMDB, taxonomy and mask files are not used.

The `.pin`, `.psq` and `.phr` files of every volume are memory-mapped with the foreign function and memory API
(`FileChannel.map` into a shared `Arena`, closed with the database); nothing is copied onto the Java heap, and the
sequence and header offsets are read in place from the mapped `.pin`. Opening a database therefore takes the same few
milliseconds whatever its size, and the total size is limited only by the address space. `SequenceDatabase` holds one
`MemorySegment` per volume (a FASTA file is one heap segment). The search copies each subject with its two sentinels
into a small buffer of its own thread (`withSentinels`), so the seeding, extension and composition code still works
on `byte[]` exactly as before.

`DatabaseLoadBench` and `BlastpSearchBench` (`source` parameter) measure the change, single thread, ms:

| | heap array (before) | mapped segments |
|---|---:|---:|
| open `blastdb/sprot` | 102 | 0.6 |
| load `uniprot_sprot.fasta` | 1,264 | 1,324 |
| HBA_HUMAN, cbs0, from FASTA / BLAST db | 2,225 / 2,068 | 2,162 / 2,072 |
| CDK1_HUMAN, cbs0 | 4,033 / 4,279 | 3,927 / 4,116 |
| BRCA1_HUMAN, cbs0 | 17,422 / 16,887 | 16,506 / 16,784 |
| HBA_HUMAN, cbs2 | 2,422 / 2,458 | 2,240 / 2,210 |
| CDK1_HUMAN, cbs2 | 4,529 / 4,446 | 4,285 / 4,334 |
| BRCA1_HUMAN, cbs2 | 17,960 / 17,614 | 17,228 / 17,310 |

Copying each subject costs nothing measurable: the differences are within the noise of the runs. The whole CLI run
of the 24 queries against the BLAST database (8 threads) takes 34 s as before, with a peak resident size of 0.72 GB
instead of 0.91 GB, most of it the mapped file in the page cache; a one-sequence query peaks at 0.28 GB instead of
0.60 GB.

The identifiers are printed as blastp prints them (`db.SeqId`, checked for every kind in `src/test/resources/blastdb`):

| Header in the FASTA file | `sseqid` | `sacc` | `saccver` |
|---|---|---|---|
| `sp\|P11111\|N1_HUMAN` | `sp\|P11111\|N1_HUMAN` | `P11111` | `P11111` |
| `ref\|NP_000509.1\|` | `ref\|NP_000509.1\|` | `NP_000509` | `NP_000509.1` |
| `gi\|123456\|ref\|NP_222222.2\|` | `gi\|123456\|ref\|NP_222222.2\|` | `NP_222222` | `NP_222222.2` |
| `pdb\|1ABC\|A` | `pdb\|1ABC\|A` | `1ABC_A` | `1ABC_A` |
| `gnl\|mydb\|xyz` | `gnl\|mydb\|xyz` | `mydb:xyz` | `mydb:xyz` |
| `lcl\|myseq` | `myseq` | `myseq` | `myseq` |
| `plainid` | `plainid` | `plainid` | `plainid` |

`stitle` is the rest of the header. A database made without `-parse_seqids` stores only the identifier
`gnl|BL_ORD_ID|n`; blastp then prints the first word of the header for all three columns, and the whole header as
`stitle`. The FASTA reader parses headers as `makeblastdb -parse_seqids` does, so a FASTA file and a database made
from it with `-parse_seqids` give the same output.

## Exhaustive baseline

`search.ExhaustiveSearch` aligns the query with every subject using cssw4j's striped Smith-Waterman. It uses BLAST's
scoring: BLAST's 11/1 gap costs are SSW's 12/1, because SSW's gap open already includes the first residue. For the
matrix, `stats.Blosum62` holds BLOSUM62 as blastp actually scores it. That matrix was measured from blastp itself: for
every pair of the 28 NCBIstdaa letters, blastp aligned `W`x10 + one letter + `W`x10 against the same with a second
letter, and the pair score is the raw score minus 220. It is not quite the NCBI `BLOSUM62` file:

* X scores -1 against every letter except `*`.
* N/B, Q/Z, W/Z and B/Z score 4, 4, -2 and 0.
* U scores like C, O like X, and J has its own row.

It reports the best alignment per subject when the Spouge E-value is at most 10.

`ExhaustiveBaseline` (in the tests) runs the 24 reference queries and compares the results with the
`-comp_based_stats 0` reference, using the best HSP per subject. Without composition-based statistics, every BLAST
HSP is a local alignment under the same scores, so the exhaustive score can only be equal or higher, and so every BLAST
hit must also be an exhaustive hit.

| | |
|---|---:|
| BLAST subjects (query-subject pairs) | 9,110 |
| found by the exhaustive search | 9,110 (all) |
| exhaustive score equal to BLAST's best HSP | 8,521 |
| ... and with the same start and end on both sequences | 8,477 |
| exhaustive score higher | 589 |
| exhaustive hits with E <= 10 | 38,322 |
| exhaustive hits missed by BLAST, counting only the first 500 per query | 408 |

* No exhaustive score is below BLAST's (the run checks this), and every BLAST subject is found.
* For 93% of the subjects BLAST finds the optimal alignment.
* Higher exhaustive scores come mostly from long multi-domain queries: FBF1 249, CFTR 179, P53 47, BRCA1 39, EGFR 34
  and MYC 30. Likely causes, not yet checked one by one: BLAST's best HSP covers only part of what Smith-Waterman
  joins into one alignment, or the X-drop extension stops early. Step 5 (gapped X-drop) has to reproduce these lower
  scores.
* 14 of the 24 queries reach BLAST's 500-subject limit. For those, only exhaustive hits ranked within the first 500
  count as missed.
* The 408 subjects BLAST misses are weak hits, close to the E-value cutoff (MYC 76, FBF1 69, BRCA1 69, P53 35, ...).
  The two-hit seeding never triggers on them, or their ungapped extension stays below the gap trigger.

These numbers are the target for the seeding steps: BLAST's sensitivity is what has to be reproduced, not the
exhaustive one.

Speed (i7-6700HQ):

Single thread, measured with `ExhaustiveSearchBench` (one search over all of Swiss-Prot):

| Query | Length | Time | GCUPS |
|---|---:|---:|---:|
| HBA_HUMAN | 142 | 7.50 s ± 0.19 | 4.0 |
| CDK1_HUMAN | 297 | 10.09 s ± 0.30 | 6.2 |

blastp on one thread takes 76 s for all 24 queries, about 3 s per query. The 24 queries add up to 10,903 residues, so
at about 6 GCUPS the exhaustive search would need about 6 minutes on one thread, roughly 5 times as long. On 8 threads
it takes 112 s (20.4 GCUPS). The 8-thread times per query vary more than their lengths explain: MYG and CRYAA took 12
and 11 s, while other short queries took 2 s. The 8-thread numbers on this laptop are also affected by thermal
throttling (see cssw4j), so they are only indicative.

**Use a HotSpot JDK (C2) for these runs, not GraalVM.** On GraalVM 25.4.4 (Oracle and CE; the sdkman default here)
a single search takes over 500 s instead of 7.5 s:

* When the Graal JIT compiles the scan loop, it leaves cssw4j's small `Lanes.shiftBytes(ByteVector)` helper out of
  line (`not-within-budget`).
* With a vector crossing that call, every vector operation in the inlined kernel then runs the Vector API's pure-Java
  fallback.
* Raising the budget (`-Djdk.graal.TypicalGraphSize=20000`) avoids it.

See the standalone reproducer in cssw4j's
[graalvm/reproducers](../cssw4j/graalvm/reproducers/README.md). All numbers above are from Temurin 25.0.4.

## Seeding

`seed.LookupTable`, `seed.WordFinder` and `search.SeedSearch` follow BLAST+'s `blast_aalookup.c`, `aa_ungapped.c` and
`blast_parameters.c`:

* **Lookup table.** For every 3-letter word, the query positions whose word scores at least 11 against it. The
  candidate words range over all 28 letters. The query's own words are always included, whatever their score.
* **Two-hit trigger.** BLAST keeps an array of the last hit per diagonal and reproduces its quirks:
  * The array length is the smallest power of two of at least query length + 40, indexed by `(q - s) & mask`.
  * A diagonal extended to the right ignores hits until they lie past the extension.
  * Instead of clearing the array for each subject, an offset is added to the stored positions.
* **Ungapped extension.** The extension goes left from the end of the second word with X-drop, and goes right only if
  it reached the first word.
* **Query-specific parameters.** BLAST computes the X-drop and the cutoff with the query's own ungapped Karlin-Altschul
  parameters. `stats.UngappedKarlin` ports that calculation: the score distribution of BLOSUM62 between the query
  composition and the Robinson and Robinson background, Newton's method for lambda, and the series for K. For all 24
  queries it gives blastp's printed lambda, K and H.
  * X-drop: `ceil(7 ln2 / lambda)`, which is 15 or 16 here.
  * Gap trigger: 22 bits, which is 40 to 42.
  * Cutoff: the gap trigger, unless the E-value cutoff gives a lower score.

**Word size 4, measured with BLAST+ only.** blastp4j seeds with 3-letter words only. A 4-letter seed only pays with a
higher neighbourhood threshold: at BLAST's 11, four letters reach the threshold so easily that the table grows and the
search gets slower. BLAST+ 2.17.0 on the 24 queries against Swiss-Prot, 8 threads, i7-6700HQ:

| `blastp` options | cycles | query-subject pairs | of the default's pairs missed | of those E < 1e-10 |
|---|---:|---:|---:|---:|
| default (word size 3, threshold 11) | 402 G | 7,817 | 0 | 0 |
| `-word_size 4` (threshold 11) | 616 G | 7,856 | 39 | 0 |
| `-word_size 4 -threshold 14` | 170 G | 7,674 | 215 | 7 |
| `-word_size 4 -threshold 16` | 81 G | 7,365 | 651 | 82 |
| `-word_size 4 -threshold 18` | 44 G | 7,023 | 1,290 | 424 |

So a 4-letter mode could only be an opt-in fast mode with different results, such as `-word_size 4 -threshold 14`
(2.4 times fewer cycles, 2.8% of the pairs missed). It would need a lookup table of its own: a word size known only at
run time made the 3-letter scan 15-27% slower (c0a3894, reverted).

`SeedingBaseline` (in the tests) runs the 24 queries and checks the result against blastp:

| | |
|---|---:|
| BLAST subjects that reach our gapped stage | 9,110 of 9,110 |
| BLAST alignments (`-comp_based_stats 0`) crossed by one of our ungapped alignments | 15,107 of 15,107 |
| blastp `-ungapped` alignments that are exactly one of ours (start, length, score) | 27,246 of 27,246 |
| Smith-Waterman hits that BLAST misses, below the 500-subject limit | 108 |
| ... of which our seeding finds no seed for | 94 |

* The second check is a requirement, not just a coincidence: a gapped extension starts at a point inside an ungapped
  alignment, so the final alignment passes through it.
* blastp `-ungapped` reports the ungapped extensions themselves. Every one of them is reproduced with the same
  coordinates and score.
* The comparison is sensitive: with the X-drop one higher, 355 of the 27,246 no longer match.
* 94 of the 108 hits BLAST misses have no seed at all. The other 14 are lost in later stages.

Speed of the seeding stage on one thread (`SeedSearchBench`, Temurin 25), against blastp on one thread:

| Query | Length | Subjects seeded | blastp4j seeding | blastp `-ungapped` | blastp (gapped) |
|---|---:|---:|---:|---:|---:|
| HBA_HUMAN | 142 | 6,339 | 1.89 s | 1.83 s | 1.96 s |
| CDK1_HUMAN | 297 | 15,883 | 3.16 s | 3.03 s | 3.65 s |
| BRCA1_HUMAN | 1,863 | 65,541 | 11.29 s | 10.52 s | 13.83 s |

The blastp times include starting the program and opening the database. Seeding takes most of blastp's time: its
whole gapped search takes only 0.1 to 3.3 s more than its ungapped one. blastp4j's seeding runs at about the speed of
BLAST's. On 8 threads, the 24 queries are seeded in 20 s.

**Query parameters as constants.** `WordFinder.of` gives a copy of `WordFinder`, made with the class-file API and
defined as a hidden class, in which the query's window, X-drop, cutoff and diagonal mask are constants (d7bdab1:
1.5-2.2% fewer cycles for CDK1 and BRCA1 in `WordFinderBench`). Single-query searches use it. The batched scan
(`BatchWordFinder`, which the CLI uses for batches of several queries) does not copy the queries' finders: all
queries' hits go through one call site, where copies would make each hit a virtual call. Since 45b95c9
`BatchWordFinder.of` copies the scan itself, with the window and the batch's diagonal mask as constants (shared by
nearly all batches: BLAST's are about 10,000 letters, so the mask is 16383); `HiddenCopies` makes both kinds of copy.
On the runners: 1.0% and 1.2% less wall time on AArch64 in two runs (every round), none on x86-64 or locally.

*Kept:* the batch scan's diagonal index also masked with the array's length (7c980c5), a power of two: C2 removes the
bounds check of `a[i & (a.length - 1)]`. 1.8% fewer instructions on the bean queries locally; on AArch64 1.7% fewer
cycles and 1.0% fewer instructions, wall time 0.3% less (no round slower); none on x86-64. *Tried, slower:* the same
in the single-query scans (498cdaa, reverted): `WordFinderBench` 1.0-1.7% slower.

*Kept:* the batch table's entries also as chars (928b663), half the memory the scan reads, where the concatenation is
shorter than 65,536 letters, as BLAST's batches of about 10,000 are. The source keeps reading ints (longer batches, the
native image); the scan's hidden copy for a narrow table is made to read chars by a further code transform
(`entry()` becomes `narrow()`, each `entry[k]` load a char load), and the verifier rejects a copy where a load was
missed, so that `of` falls back to the ints. `BatchWordFinderTest` checks both copies against the int scan. On the
runners: 0.8% less wall time on AArch64 and 1.4% on an EPYC 9V45, in every round; locally within the noise. (Chars
for a single query's table, 586aa09, gained nothing: that table fits in L1.)

*Tried, no gain:* the single-query copies with the batch's diagonal base (0) and its cross-query check (never true for
a query alone) as constants too (48bc96d, reverted): `WordFinderBench` unchanged, as the check is on the rare path.

*Tried, slower:* a copy of `BatchWordFinder` per thread that sends each hit through a switch on the query number to
that query's copy, held as a constant (class data) (bb46e1d, reverted). C2 does not inline the 24 cases into the scan,
so every hit became a call to a separately compiled `probe`: 37% more cycles in the full search.

*Kept:* the merged table's entries as one int each, query number and offset packed, as BLAST's table for a batch holds
them (7eedf0a). Before, the query and the offset were in two arrays, so a hit could touch two cache lines: 3.6% fewer
cycles in the full search (median of eight interleaved rounds). *Not built:* BLAST's backbone, with up to three entries
inline in each word's cell. In the merged table of the 24 queries a non-empty word has 22.7 entries on average, and
only 493 of 15,666 have three or fewer, so the inline entries would rarely be used.

*Kept:* one diagonal array for the batch, as BLAST has for its concatenated queries. Each query keeps its own
power-of-two number of diagonals, in a region that starts at a multiple of its size. An entry holds that size's
logarithm and the region start plus the query offset, so the batched scan finds a hit's diagonal from the entry
alone, and handles the usual hit (the diagonal's last hit a window or more back) without touching the query's
`WordFinder`. The query number is in a parallel array, read only for the other hits. Full search, three interleaved
rounds: 2.1% fewer cycles and 3.4% fewer instructions on the 24 human queries against Swiss-Prot, 1.4% fewer cycles
on the 500 bean queries.
On the GitHub runners (2755a0e): 2.9% less wall time and cycles on AArch64, 1.6% less wall time on x86-64.
Since b9735a8 the array is exactly BLAST's: indexed with the offset in the concatenated queries (each query one
separator after the previous one, as `s_QueryInfo_SetContext` places them) and sized for the whole concatenation, so
an entry is that offset and its diagonal `(entry - s) & mask`, without decoding a region. Diagonals then run across
queries, and a hit whose diagonal's last hit lies in the previous query starts a new pair, as in BLAST (its check
after `BSearchContextInfo`); without that check one line of `BlastpTest.optionsAsBlastp` differed from blastp. Full
search, local, against the regions: 2.9% fewer cycles on the bean queries, 4.6% on the human ones (5.9% and 8.2%
fewer instructions).

*Kept:* the usual hits decided without branching on them (493154b batched, 59b0741 single query). On the bean
queries a quarter of the 5.0 G entries the batched scan visits were not the usual store; three quarters of those only
overlapped the diagonal's last hit and were ignored, yet each went through the query's `WordFinder` first. Now the scan
ignores them itself, and one sign test decides between storing, ignoring and calling the finder, with a masked store.
Written with `||`, C2 compiled the test to separate branches and nothing changed; with the sign test, branch misses
fell by 24% and cycles by 5.4% (bean) and 4.1% (human), and `WordFinderBench` took 9-10% less time. On the runners:
6.0% less wall time on x86-64 (a Xeon 6973P-C, now 7.5% faster than BLAST+) and 3.0% on AArch64, where branch misses
fell by 20% but cycles by only 0.8%.

*Tried, no gain:* the hits for the finder put in a buffer and handled after each position's entries (7a12106,
reverted), so that the entry loop has no call. C2 still kept the loop's values in stack slots: 2.8% more instructions,
cycles unchanged.

*AArch64:* C2 has no callee-saved registers there, so the values a loop keeps across a call live in stack slots, and
the batch scan reloaded them after every call to the finder; GCC keeps BLAST+'s in x19 to x28. On the Neoverse-N2
runner 14.3% of the samples were on stack loads in `find`. Experiments, 500 bean queries, output identical:

| Commit | Change | AArch64 (Neoverse-N2) | x86-64 | |
|---|---|---:|---:|---|
| 600a57a, e66bdfe | 7a12106's deferral again: a position's rare hits buffered and handed to the finder after its loop | -3.7% wall, stack-load samples 12.0% | +2.8% (Xeon 6973P-C) | kept on AArch64 only (`BatchWordFinder.DEFER`) |
| c2ee363 | `WordFinder.hit` inlined into the scan (`jvm.options`, AArch64) | +0.3% wall, +3.5% cycles | | reverted (ac8b292) |
| c28fe14 | `seed.SeedList`: the alignments in an array grown by 64, so that saving one is no call | +0.2% wall, +1.0% cycles | -0.2% (EPYC 9V74) | reverted (f4f6445) |
| cfdae36 | 654b5e5 again: the ungapped extensions' best score without a branch | -1.2% wall, -5.7% cycles, -8.5% branch misses | (+1.8% in 654b5e5) | kept on AArch64 only (`WordFinder.BRANCH_FREE_MAX`) |
| 7c4ce9c | a position loop with no call at all (`scan`), returning to `find` for the hits | +0.1% wall, +2.0% cycles | | reverted (0d692cf) |

Without the inline and the array (their reverts together) the N2 took 4.4% less wall time, every round. The call-free
scan did remove the reloads (0.2% of the samples on stack loads), but stored a value per position and made `find`
6,280 instructions: as a8ca80b found, the samples on the reloads were mostly the latency of the table and diagonal
loads before them. Both constants are `os.arch` tests, folded by C2; `seed.FastWordFinder` has the branch-free
maximum too.

*Kept:* the vector scan of a single query (`WordFinder.findBlocks`) at most 256 bits wide (1aff999). With AVX-512 it
took 16 positions per vector: on a Xeon 8573C a search of CDK1 took 10% longer with them than with 8
(`BlastpSearchBench`, `-XX:MaxVectorSize=32` against the default, two runs), and on an EPYC 9V45 (Zen 5) the full
search of the 500 bean queries took 7.6% less wall time with the cap, which brought it to BLAST+'s (33.1 s, BLAST+
33.6 s). Machines without AVX-512 already used 8 lanes or fewer. *Since:* the code uses the preferred width everywhere
and `bin/jvm.options` chooses it: 256 bits on x86-64 (`-XX:MaxVectorSize=32`, also on Zen 4 and Intel's AVX-512
Xeons), 512 bits on AMD Zen 5 (`[zen5]`: `-XX:UseAVX=3 -XX:MaxVectorSize=64`), there for the word scan too (ea24369,
d17f101). Measured since on an EPYC 9V45 (main's push run 42221716347, against 752525a): 29.33 s wall against
30.15 s (-2.7%), BLAST+ 33.82 s, almostblastp 27.89 s; output identical. The earlier loss with 16 lanes did not come
back: the scan has changed much since.

## Gapped extension and traceback

`align.GappedAligner`, `align.HspTree` and `search.BlastpSearch` follow BLAST+'s `blast_gapalign.c`, `blast_itree.c`,
`blast_hits.c`, `blast_traceback.c` and `blast_engine.c`. Per subject:

1. **Preliminary gapped extension.** For each ungapped alignment, best first, that is not inside a gapped one already
   found:
   * The start point is the middle of the alignment's best-scoring 11-letter window.
   * The extension is score only, with X-drop 15 bits (38 raw).
   * If the subject's best ungapped score is below 0.68 times the cutoff, the restricted, approximate alignment is used
     (gaps may only start every 10 positions). When its result is inconclusive, everything for that subject is redone
     exactly.
2. **Preliminary hit list.** HSPs sharing a start or end point are purged, E-values are computed, and the best 550
   subjects are kept.
3. **Traceback.** Each HSP not inside an already traced one is aligned again from its start point, with traceback and
   X-drop 25 bits (64 raw). Leading and trailing gaps are trimmed, then purging, containment and the E-value cutoff
   follow. The final list is the best 500 subjects.

`GappedBaseline` (in the tests) runs the 24 queries and compares every output line with the `-comp_based_stats 0`
reference:

| | Lines identical to blastp's (of 15,107) |
|---|---:|
| first complete version | 14,926 |
| + BLAST's interval-tree endpoint rule | 15,093 |
| + BLAST's minimum subject length | **15,107**, in the same order |

Three details took the result from 98.8% to exact. All three are visible only in BLAST's source or behaviour:

* **The containment tree is not just a set.** Before a gapped alignment is computed, BLAST checks whether the ungapped
  alignment lies inside an HSP already found. When an HSP is added to that tree, any HSP already there with the same
  start (or end) point is removed if it scores lower. If the one already there is better, the new one is not added.
  So a lower-scoring HSP that was replaced no longer hides the alignments inside it. `HspTree` reproduces this.
* **The minimum subject length is 10.** The score an HSP needs to be kept is derived from the E-value cutoff at the
  shortest subject length. BLAST databases do not store that length, so SeqDB reports 0 and BLAST uses its default
  `BLAST_SEQSRC_MINLENGTH` of 10. This cutoff also decides when the restricted alignment is used, so the real
  minimum (2 in Swiss-Prot) changed weak hits near the cutoff.
* **Numbers are rounded as C's printf does.** printf rounds the exact binary value with ties to even, so 25/64 =
  39.0625 prints as `39.062`. Java's `String.format` rounds ties up.

Some BLAST quirks are kept on purpose, because the output depends on them:

* When both extensions end in a gap of the same kind at the start point, the edit script merges the two gaps into one,
  but the score still counts both gap openings.
* An HSP whose traceback start is exactly (0, 0) gets a new start point.
* The score-only and the traceback extension disagree by one letter about the end of an alignment that starts at the
  last letter of a sequence.

`BlastpSearchReferenceTest` checks three of the queries line by line in `mvn test`, when Swiss-Prot is present.
`GappedAlignerTest` checks on random sequences that the traceback agrees with the score-only extension and with its
own edit script, and that the restricted alignment never scores higher than the exact one.

Speed of the complete search on one thread (`BlastpSearchBench`, Temurin 25), against blastp on one thread:

| Query | Length | blastp4j | blastp | blastp4j seeding only |
|---|---:|---:|---:|---:|
| HBA_HUMAN | 142 | 2.10 s | 1.96 s | 1.89 s |
| CDK1_HUMAN | 297 | 4.13 s | 3.65 s | 3.16 s |
| BRCA1_HUMAN | 1,863 | 17.28 s | 13.83 s | 11.29 s |

The blastp times include starting the program and opening the database. blastp4j is 7 to 25% slower. The difference
is in the gapped stage: for BRCA1 it takes 6.0 s here, against about 3.3 s in blastp (its gapped search minus its
ungapped one). On 8 threads, the 24 queries take 29 s; blastp on one thread takes 76 s.

*Tried, slower:* the score-only extensions in 4 lanes on AArch64 with SVE (80b15ff, reverted), where they run scalar
because the lanes need 8 ints with a gather. On the Neoverse-N2 runner 6.2% more wall time than scalar, though with 4.6%
fewer instructions; on x86 4 lanes take 15% more cycles than scalar, 8 lanes 1% fewer.

*Kept:* the traceback's row loop as a hidden-class copy per pair of gap costs (`GappedAligner.ExRows`, dfc3c67), as the
score-only loop has. Before, the inlined loop had constant gap costs only for (11, 1); with composition-based
statistics, the default, the traceback aligns with everything scaled by 32, so it ran with the costs as variables (the
runners' machine-code listings showed them in stack slots). On the runners: 0.5% less wall time on AArch64, with 1.3%
fewer cycles and 2.2% fewer instructions; none on x86-64. A third inlined call for the scaled costs instead (dada6ae,
reverted) gained nothing. `RowsTest` compares the copies' tracebacks with `ExRows` itself, also with scaled costs.

*Tried again, no gain* on top of the changes since (one run each on the runners, 500 bean queries, 3 rounds):

| Retry of | Change | AArch64 | x86-64 |
|---|---|---|---|
| fa4154f | traceback's best and bestGap interleaved (fa20928) | -0.1% | -0.1% |
| eee7c94 | traceback's X-drop test against a floor (24f7c19) | rounds mixed | (job failed) |
| c86b5b3 | ungapped extension without a branch on a new maximum (654b5e5) | -0.8% | +1.8% |
| e9d5cb1 | ungapped extension loops unrolled by four (54574b8) | +2.3% | +4.7% |
| 36f465f | score-only row's last live cell found after the row (aa1ef5b) | +2.4% | rounds mixed |

*Kept:* the composition-based stage's tracebacks in lanes, many subjects at once (`align.TracebackLanes`,
`compo.TracebackStage`; branch `exp/traceback-lanes`, ea8ac72 to c08e2df), where `SemiGappedLanes` runs (x86-64,
preferred vectors of 8 ints or more). Each lane runs one side of one traceback, alignEx's dynamic program with a
traceback byte per cell and the subject's own adjusted matrix; 16 subjects per task are in flight, each driven by
`CompoRedo.Subject`, the scalar redo's steps. Alone the kernel took 1432 ms where the scalar one took 2848 ms. On the
runners (500 bean queries): x64 wall -10.2% on an EPYC 9V74 and -8.2% on an EPYC 7763, AArch64 unchanged (the lanes
are off there; the scalar redo is main's own loop, as driving `Subject` there cost 7% more instructions). Locally 7.8%
fewer cycles; memory up about 15% (the trace rings).

*Tried, slower:* on AArch64, where the lanes are off, two scalar tracebacks interleaved cell by cell
(`align.TracebackPairs`, 7d151d5, reverted in 66eeae5), so that the core has a second dependency chain: alignEx's
cells are one serial chain, gap row to score to gap row. On the Neoverse-N2 (run 42469313465) blastp4j took 7% more
wall time and fastblastp 6% more: the pair kernel had 139k samples where `ExRows` had 102k, its fused loop full of stack
loads (two extensions' row state does not fit in AArch64's registers under C2), and it has none of `ExRows`'
score-only rows after the best.

*Measured, not built:* score-only extensions in 16 lanes of shorts instead of 8 ints, for twice the cells per vector
(an idea for almostblastp after the G-BLASTN paper). There is no 16-bit gather, and the score matrix is read with a
gather because both letters differ per lane, so 16 lanes need two 8-int gathers and a narrowing to shorts per slot;
only the arithmetic and the DP row's memory halve. The bound, `SemiGappedLanes.row`'s cell loop alone with
saturating adds (a scratch kernel, i7-6700HQ, AVX2, JDK 27 C2, warm): 1.08 to 1.12 ns per cell with 16 short lanes
against 0.94 to 0.98 ns with 8 int lanes, so 15% slower per cell before any bookkeeping. Not pursued.

*Tried, kept as branches only:* ungapped extensions in lanes across hits (`exp/ungapped-lanes`, c2aecec, 1200971):
40-85% more cycles locally (plain and vector lanes), 45-52% more wall time on the runners. The traceback filter (`exp/traceback-filter`, b881f5c): a proven bound on the traceback score from a
score-only extension, so that a subject with one preliminary alignment whose bound cannot reach the E-value cutoff
skips its traceback; -1.9% cycles locally, -0.3% to -0.4% wall on the runners, within noise.
The traceback row loop's stack reloads on AArch64 (`ExRows.rows` on the Neoverse-N2 reloads 7 values from the stack
per cell): the best score, its X-drop floor and column in locals (`exp/exrows-registers`, c02ad8f), then `best` and
`bestGap` interleaved in one array (`exp/exrows-interleaved`, 9a1f896), then the matrix row copied per row
(`exp/exrows-rowscores`, 8e5f20b). Against b6e51a2 on the ARM runner: 6, 5 and 5 reloads per cell; instructions
-1.3%, -2.5% and +3.9%; cycles +0.5%, 0.0% and +2.8%; wall within noise on both runners. Fewer reloads save
instructions but not cycles: the reloads are not what limits this loop.
*Kept:* the cell loop in a method of its own (`ExRows.row`, dfb86fb), kept out of `rows` with a `CompileCommand
dontinline`, so that C2 allocates its registers for the loop alone. The hottest reload went from 4.3% to 1.9% of the
fastblastp samples on the N2, and `ExRows` from 18.3% to 17.0%. Wall time (2,500 queries, median of three rounds),
run 42510201537 and its rerun 42513867132: on the N2, blastp4j +0.2% and -0.7%, fastblastp -3.2% (noisy reference)
and -0.5%; on x86-64, blastp4j -0.9% (EPYC 9V74) and +0.4% (EPYC 9V45), fastblastp +0.2% and +0.5%. Locally
(i7-6700HQ, scalar path) instructions +1.6%, cycles neutral.
After it, `ExRows.AHEAD` (rows with traceback after the best cell's before the score-only rows) again: 8 rows
instead of 4 (a9c6886) is kept, fastblastp -2.1% and -2.4% on the N2 (all six rounds faster), blastp4j +0.6% and
-0.4%, x86-64 within 0.4% (runs 42518240547, 42521981040). 2 rows (cac2447) is reverted (8290c0b): on the N2
fastblastp +3.9% and +1.5%, blastp4j 0% and +1.2%; x86-64's blastp4j -5.1% in the first run was not there in the
rerun (-0.3%).

## Composition-based statistics

blastp's default, mode 2 ("conditional compositional score matrix adjustment", Altschul et al. 2005; Yu and
Altschul 2005), changes the search in three places:

* **Preliminary stage.** The cutoffs use 5 times the E-value cutoff (50), and 1,050 subjects are kept instead of 550.
* **Traceback stage.** It is replaced. For each subject (`Blast_RedoAlignmentCore` and `Blast_RedoOneMatch` in
  `blast_kappa.c` and `redo_alignment.c`):
  * The subject is masked with SEG (window 10, locut 1.8, hicut 2.1), unless its first alignment shows it to be nearly
    identical to the query.
  * The amino-acid compositions of query and subject decide how BLOSUM62 is adjusted. Usually its target frequencies
    are re-optimized for the two compositions at a relative entropy of 0.44 (Newton's method on a constrained
    Kullback-Leibler problem). When the compositions differ a lot and the lengths differ more than threefold, the
    matrix is instead only rescaled ("old" composition-based statistics).
  * The preliminary alignments are aligned again with that matrix. All scores are scaled by 32 for precision.
* **Result list.** BLAST reads the subjects in database order and keeps them in a heap that holds every subject
  below E = 0.002 and otherwise the best 500. Once the heap is full of such subjects, a subject whose preliminary
  E-value is above 0.01 is skipped.

The port, in `swiss.sib.blastp4j.compo`:

* `Seg` is BLAST's SEG with its tabulated logarithms. It masks exactly what BLAST's `segmasker` does, on the first
  20,000 Swiss-Prot sequences (`SegTest` uses a sample of 400).
* `CompositionAdjustment` and `OptimizeTargetFreq` compute the adjusted matrices. `Blosum62Data` holds BLAST's
  BLOSUM62 joint probabilities and frequency ratios.
  * All matrices of the stage are one array each, row-major (b86b14c; the Newton matrix of `OptimizeTargetFreq` a lower
    triangle packed by rows), for the simpler layout, not for speed: `CompoRedoBench` was 1.3% slower. `ScoreMatrix`
    holds a matrix for the aligner, immutable, so a subject's adjusted matrix is handed on without copying.
  * `bin/jvm.options` keeps C2 from inlining `OptimizeTargetFreq`'s methods (a2a0da3). With the flat arrays, inlining
    was tried again (5b96a04, reverted): 1.3% (AArch64) and 1.6% (x86-64) more wall time on the GitHub runners, 4.8%
    more cycles here.
  * The Newton system's matrix A D A^T has a diagonal block, the first 20 rows (the column constraints). BLAST's
    Cholesky factorization and solve run over its zeros as over any entry; `OptimizeTargetFreq` leaves them out
    (8dff3b7), about half the factorization's multiply-adds. Subtracting the products with them changes nothing, and
    the other terms keep BLAST's order, so the results are BLAST's. On the runners: 1.7% less wall time on AArch64
    (2.9% fewer cycles) and 0.5% on an EPYC 9V45, every round; 2.3% fewer cycles on the bean queries here.
  * *Tried, no gain:* the factorization and solve with fused multiply-adds, four partial sums and the diagonal's
    reciprocals, for almostblastp (39f9f3e, reverted): the same results, but no faster; the dot products are short.
  * BLAST's own C code was compiled into a test harness: `reference/compo/build_harness.sh` downloads the library and
    builds `harness.c`.
  * Of 6,006 query/subject compositions, every adjusted matrix is identical to the C code's. One λ ratio differs in
    its last bit, which no output depends on. `CompositionAdjustmentTest` checks 306 of them.
* `CompoRedo` redoes the alignments of one subject: the near-identity tests, SEG, the adjustment, the realignment,
  BLAST's rules for distinct and contained alignments, and the scaled E-values.
* `CompoHeap` is BLAST's result heap.

The realignment runs in parallel. Early termination and the heap are then replayed in BLAST's order, which gives the
same result because one subject's realignment does not depend on another's. There is one exception, and it is
reproduced. When SEG leaves a subject with no true amino acids, the adjustment fails, and BLAST realigns that
subject's later alignments with the matrix of the previous subject it processed. The replay then redoes such a
subject with that matrix.

Two details were needed for an identical order of the output lines:

* E-values are computed as BLAST does, from the scaled score with λ / 32 and the Gumbel parameters a, α and σ scaled
  by the same ratio. Dividing the score by 32 gives the same value mathematically, but not to the last bit.
* BLAST's C++ layer sorts the HSPs of each subject by E-value when it converts the results to alignments
  (`blast_seqalign.cpp`). Because E-values below 1e-180 count as equal, HSPs with the same final score are then
  ordered by position, not by their scaled score.

`GappedBaseline` with the argument `cbs2` compares all 24 queries: 11,263 of 11,263 lines are identical, in the same
order. `BlastpSearchReferenceTest` checks three of them in both modes.

Speed of the complete search on one thread (`BlastpSearchBench`, Temurin 25), against blastp on one thread:

| Query | blastp4j, mode 2 | blastp, default | blastp4j, cbs0 | blastp, cbs0 |
|---|---:|---:|---:|---:|
| HBA_HUMAN | 2.43 s | 2.23 s | 2.10 s | 1.96 s |
| CDK1_HUMAN | 4.45 s | 3.89 s | 4.13 s | 3.65 s |
| BRCA1_HUMAN | 18.16 s | 14.07 s | 17.28 s | 13.83 s |

On 8 threads, the 24 queries take 32 s with composition-based statistics; blastp on one thread takes 80 s.

## Reference results

See [reference/README.md](reference/README.md): BLAST+ 2.17.0 on UniProtKB/Swiss-Prot 2026_03 for 24 human proteins,
with and without composition-based statistics, regenerated by `reference/make_reference.sh`.

The reference is BLAST+ on Linux. BLAST+'s own output is not the same everywhere
([docs/blast-determinism.md](docs/blast-determinism.md)):

* **Platform.** On macOS, BLAST+ differs on 2 of 500 bean queries. Floating point is not the cause; the order of ties
  in `qsort` is. glibc's sort is stable and Apple's is not. Preloading a sort with Apple's order
  (`tools/qsort-order/apple-qsort.c`) into Linux BLAST+ reproduces the macOS output byte for byte. blastp4j sorts
  stably, so it gives Linux's output on macOS too.
* **Threads.** BLAST+'s output does not change with the thread count unless a query fills its 500-subject heap while
  alone in its batch. blastp4j gives the single-threaded result.
* **Batches.** A different `BATCH_SIZE`, or a different set of queries in a batch, can change other queries' gapped
  results (10 of 500 queries with all of them in one batch). The likely cause is BLAST's interval tree over the
  concatenated batch. blastp4j matches the default batches only.

## Build

Requires JDK 25 or newer and cssw4j installed in the local Maven repository (`mvn install` in the cssw4j project). The
classes are compiled for Java 25, so that GraalVM 25's native-image can build the native blastp4j. Runs and benchmarks
use JDK 27 (C2) as the reference, as does the architecture comparison workflow; numbers measured on Temurin 25 say so.
That workflow (`.github/workflows/arch-compare.yml`) runs BLAST+ and blastp4j on x86-64 and AArch64 runners and
publishes a report of each run, with each runner's CPU model, at <https://www.jerven.eu/blastp4j/>.

    mvn test

The baselines and their benchmarks need Swiss-Prot in `~/data/uniprot` (see `reference/make_reference.sh`):

    mvn test-compile dependency:build-classpath -Dmdep.outputFile=target/cp.txt
    CP=target/classes:target/test-classes:$(cat target/cp.txt)
    java --add-modules jdk.incubator.vector -Xmx3g -cp $CP swiss.sib.blastp4j.search.ExhaustiveBaseline 8
    java --add-modules jdk.incubator.vector -Xmx3g -cp $CP swiss.sib.blastp4j.search.SeedingBaseline 8
    java --add-modules jdk.incubator.vector -Xmx3g -cp $CP swiss.sib.blastp4j.search.GappedBaseline 8 cbs0
    java --add-modules jdk.incubator.vector -Xmx3g -cp $CP swiss.sib.blastp4j.search.GappedBaseline 8 cbs2
    java -cp $CP org.openjdk.jmh.Main 'ExhaustiveSearchBench|SeedSearchBench|BlastpSearchBench'

`SeedingBaseline` also counts the Smith-Waterman hits without a seed when `target/exhaustive.tsv` from
`ExhaustiveBaseline` is present.
