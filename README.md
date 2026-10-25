# blastp4j

A protein database search in Java that aims to reproduce the output of NCBI BLAST+ `blastp`: the same hits, raw
scores, bit scores, E-values and alignments, for the same parameters. It builds on
[cssw4j](../cssw4j) (a Java Vector API port of the SSW Smith-Waterman library) and takes ideas from
[BLVector](https://pmc.ncbi.nlm.nih.gov/articles/PMC7884812/).

Status: the statistics (bit scores and E-values, reproduced exactly for all 15,107 hits of the reference without
composition-based statistics) and the database and query reading (identical to BLAST's own database, byte for byte)
are implemented, as is an exhaustive Smith-Waterman baseline search. The BLAST heuristic search is not implemented yet.

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
   1.1 s from plain FASTA and 2.3 s from gzip. Reading BLAST database files directly is left for later.
3. **Exhaustive baseline** (done) with cssw4j over all of Swiss-Prot; see [the results](#exhaustive-baseline).
4. **Seeding**: lookup table, two-hit trigger, ungapped X-drop extension. Check: which subjects reach the gapped stage.
5. **Gapped X-drop extension and traceback.** Check: raw scores, coordinates and `btop` of the
   `-comp_based_stats 0` reference.
6. **Composition-based statistics** (mode 2). Check: the default (`cbs2`) reference.
7. **Command line, tabular output (`-outfmt 6`), threads and hit list rules**; speed against blastp.

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

**Use a HotSpot JDK (C2) for these runs, not GraalVM.** On Oracle GraalVM 25.4.4 (the sdkman default here), the
Graal JIT recompiles cssw4j's 8-bit kernel after a deoptimization. In the recompiled code the Vector API calls are not
intrinsified: the stack shows the pure-Java fallback `ByteVector.bOpTemplate`. A single search then takes over 500 s
instead of 7.5 s. All numbers above are from Temurin 25.0.4.

## Reference results

See [reference/README.md](reference/README.md): BLAST+ 2.17.0 on UniProtKB/Swiss-Prot 2026_03 for 24 human proteins,
with and without composition-based statistics, regenerated by `reference/make_reference.sh`.

## Build

Requires JDK 25 and cssw4j installed in the local Maven repository (`mvn install` in the cssw4j project).

    mvn test

The step 4 baseline and its benchmark need Swiss-Prot in `~/data/uniprot` (see `reference/make_reference.sh`):

    mvn test-compile dependency:build-classpath -Dmdep.outputFile=target/cp.txt
    CP=target/classes:target/test-classes:$(cat target/cp.txt)
    java --add-modules jdk.incubator.vector -Xmx3g -cp $CP swiss.sib.blastp4j.search.ExhaustiveBaseline 8
    java -cp $CP org.openjdk.jmh.Main ExhaustiveSearchBench
