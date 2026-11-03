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
(`BatchWordFinder`, which the CLI uses for batches of several queries) does not: all queries' hits go through one call
site, where copies would make each hit a virtual call.

*Tried, slower:* a copy of `BatchWordFinder` per thread that sends each hit through a switch on the query number to
that query's copy, held as a constant (class data) (bb46e1d, reverted). C2 does not inline the 24 cases into the scan,
so every hit became a call to a separately compiled `probe`: 37% more cycles in the full search.

*Kept:* the merged table's entries as one int each, query number and offset packed, as BLAST's table for a batch holds
them (7eedf0a). Before, the query and the offset were in two arrays, so a hit could touch two cache lines: 3.6% fewer
cycles in the full search (median of eight interleaved rounds). *Not built:* BLAST's backbone, with up to three entries
inline in each word's cell. In the merged table of the 24 queries a non-empty word has 22.7 entries on average, and
only 493 of 15,666 have three or fewer, so the inline entries would rarely be used.

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
