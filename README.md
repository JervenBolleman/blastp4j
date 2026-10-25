# blastp4j

A protein database search in Java that aims to reproduce the output of NCBI BLAST+ `blastp`: the same hits, raw
scores, bit scores, E-values and alignments, for the same parameters. It builds on
[cssw4j](../cssw4j) (a Java Vector API port of the SSW Smith-Waterman library) and takes ideas from
[BLVector](https://pmc.ncbi.nlm.nih.gov/articles/PMC7884812/).

Status: the statistics (bit scores and E-values, reproduced exactly for all 15,107 hits of the reference without
composition-based statistics) and the database and query reading (identical to BLAST's own database, byte for byte)
are implemented. The search itself is not implemented yet.

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
3. **Exhaustive baseline** with cssw4j over all of Swiss-Prot: which reference hits it finds, how its scores compare,
   and how fast it is.
4. **Seeding**: lookup table, two-hit trigger, ungapped X-drop extension. Check: which subjects reach the gapped stage.
5. **Gapped X-drop extension and traceback.** Check: raw scores, coordinates and `btop` of the
   `-comp_based_stats 0` reference.
6. **Composition-based statistics** (mode 2). Check: the default (`cbs2`) reference.
7. **Command line, tabular output (`-outfmt 6`), threads and hit list rules**; speed against blastp.

## Reference results

See [reference/README.md](reference/README.md): BLAST+ 2.17.0 on UniProtKB/Swiss-Prot 2026_03 for 24 human proteins,
with and without composition-based statistics, regenerated by `reference/make_reference.sh`.

## Build

Requires JDK 25 and cssw4j installed in the local Maven repository (`mvn install` in the cssw4j project).

    mvn test
