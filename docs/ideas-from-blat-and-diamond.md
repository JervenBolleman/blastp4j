# Ideas from BLAT and DIAMOND

Notes on two fast protein aligners and what blastp4j might take from them. blastp4j's own constraint is that
`blastp4j` gives BLAST+'s output exactly, while `almostblastp` may differ slightly (`README.md`). So each idea is marked
**exact** (can stay BLAST's output, a candidate for both programs), **almost** (changes results a little:
almostblastp only) or **different program** (changes which hits are found: neither, as things stand).

Sources:

* W. J. Kent, *BLAT — The BLAST-Like Alignment Tool*. Genome Research 12:656–664 (2002).
* B. Buchfink, K. Reuter, H.-G. Drost, *Sensitive protein alignments at tree-of-life scale using DIAMOND*. Nature
  Methods 18:366–368 (2021), doi:10.1038/s41592-021-01101-x (v2.0.7; full text via Europe PMC, PMC8026399).

Where blastp4j's time goes now (500 bean queries against the bean proteome, AArch64 runner, share of samples): the
batch word scan about 25%, the traceback row loop 19%, the score-only gapped extension 15%, the ungapped extensions
14%, the composition-based stage's optimizer about 5% (almostblastp: vectorized), SEG now cached.

## BLAT

BLAT indexes the **database** (non-overlapping K-mers, built once, kept in memory) and scans the **query**, the reverse
of BLAST. For proteins it seeds with one perfect 5-mer by default (three perfect 4-mers in the client/server
version), and analyses at length how sensitivity and the number of chance hits trade off between single perfect,
near-perfect (one mismatch) and multiple perfect K-mer matches near each other on a diagonal (Tables 3–10: for
proteins, two perfect 4-mers within a window give BLAST-like sensitivity with far fewer chance hits than one perfect
5-mer). Hits are collected, put in buckets of 64 k database positions and **sorted on the diagonal** (database minus
query position), then clumped; protein alignment chains ungapped HSPs through a graph instead of a gapped dynamic
program.

* **Index the database, scan the query** (*exact*, a large redesign). blastp4j scans the whole database for every
  query batch. For the daemon, which serves many jobs against the same database, a database index would turn each
  search into lookups of the query's words. To keep BLAST's output, the hits of each subject must then be processed
  in BLAST's order (by subject position, then query offset), so they would be sorted per subject before the
  two-hit logic: the cost to measure against the scan. Index size is the catch: one entry per database position
  (Swiss-Prot 209 M letters, so 1–2 GB with 4–8 bytes per entry) and BLAST's neighbourhood words (threshold 11)
  must be expanded on the query side, which BLAST does already.
* **Hits sorted on the diagonal before the two-hit logic** (*exact*, a candidate for the batch scan). See the same
  idea in DIAMOND below.
* **Seed criteria** (two perfect 4-mers, near-perfect 8-mers, one perfect 5-mer) (*different program*): BLAST's
  neighbourhood words of length 3 with the two-hit rule are what blastp4j must reproduce.
* **Chaining ungapped HSPs instead of a gapped extension** (*different program*): it loses indels (Kent notes it).

## DIAMOND (v2)

DIAMOND's search stage is **double indexing**: seed–location tables are built on the fly for both the queries and a
chunk of the reference, and matched with a hash join (recursive radix clustering until a cluster's query table fits
in cache), "avoiding the memory latency bottleneck of a classical single-indexed seed lookup"; with up to 64 spaced
seed shapes (weight 7 to 10), one at a time, the indices and the join take under 1% of the run time. After the join
comes a chain of vectorised filters, each passing on fewer hits: a Hamming-distance filter on 48-letter windows
(SSE `pcmpeqb`/`pmovmskb`/`popcnt`, ten times faster than an ungapped extension with a matrix), ungapped extensions
of one query against **32 subjects at once** (AVX2, after a 32 × 32 byte transposition), a leftmost-seed filter, an
**adaptive ranking** of targets by ungapped score (extensions in chunks of 400 targets, stopping when a chunk yields
nothing significant), a **gapped extension filter** (ungapped diagonal scores combined by a one-dimensional dynamic
program, an overestimate that cuts spurious hits tenfold), **chaining** of diagonal segments to place a band, and
finally a **banded SWIPE**: one query against up to 32 subjects, 8-bit scores, recomputed on overflow. It can apply
BLAST's composition-based matrix adjustment (`--comp-based-stats`), masks with tantan instead of SEG by default, and
reaches BLASTP's sensitivity in `--ultra-sensitive` mode at 80–360 times its speed.

* **Cache-aware seed matching: sort or radix-cluster the hits, then process them in order** (*exact*, the most
  relevant idea for blastp4j's largest cost). The batch scan's time per table entry is dominated by two random loads,
  the word's cell and the diagonal's state (`README.md`, the machine-code listings). Per subject, the scan could
  instead emit (diagonal, subject position, entry) triples, radix-sort them by diagonal (in fact by BLAST's diagonal
  array index, `(entry - s) & mask`, so that diagonals sharing an entry still meet in order) and then by position,
  and run the two-hit logic one diagonal at a time with its state in a register. Each diagonal's hits are then
  handled in BLAST's order; the alignments found per query must be put back in BLAST's order afterwards (by the
  second hit's position, then query offset). Whether the sort costs less than the random loads it removes (about 21
  entries per position on the bean queries) is the experiment.
* **Ungapped extensions across hits, many at once** (*exact*). DIAMOND vectorises the ungapped extension across 32
  subjects, not along one extension. blastp4j's ungapped extensions (`WordFinder.extendLeft` and `extendRight`,
  14%) run one hit at a time; vectorising along one extension was slower (README: bb9579d), but lanes of different
  hits, each with its own X-drop, keep each extension BLAST's (integer arithmetic). The difficulty is the order: an
  extension's result decides the diagonal's state for later hits, so hits would be extended in batches only where
  they are independent (different diagonals), with the state updates applied in order afterwards.
* **The traceback across alignments, many at once** (*exact*). DIAMOND's final extension is SWIPE, inter-sequence
  SIMD. blastp4j already does this for the score-only extension on x86-64 (`SemiGappedLanes`, 8 lanes), not for the
  traceback (`GappedAligner.ExRows`, 19–22% of the time), which in the composition-based stage realigns each kept
  alignment of a subject with its adjusted matrix. Lanes would need a traceback byte per cell and lane, and the
  realignments of one subject depend on each other (containment, distinct ends), so lanes would hold different
  subjects. On AArch64 4-lane score-only extensions were slower than scalar (README: 80b15ff), so the gain there is
  doubtful.
* **Adaptive ranking and the gapped extension filter** (*almost* at best, more likely *different program*): both
  skip extensions that BLAST performs, so a hit BLAST reports could be lost. Related, and possibly exact: the
  composition-based stage spends about half its tracebacks on subjects whose result is empty (measured: 79,746 of
  164,781 on the bean queries, but only 12% of the aligned letters); skipping them needs a proven upper bound of the
  traceback score, which this filter (an overestimate) could provide. Queued in README as the traceback filter.
* **8-bit scores with recomputation on overflow** (*exact*): narrower lanes double the lanes per vector. blastp4j's
  scaled scores (composition-based, × 32) do not fit 8 or 16 bits; the unscaled preliminary stage's might in 16 bits.
* **Spaced seeds, tantan masking, chaining, banded extension** (*different program*).

## Order of experiments

1. Hits sorted by diagonal in the batch scan (exact; the largest cost).
2. The composition stage's traceback filter with a proven upper bound (exact; up to about 12% of the traceback work).
3. Ungapped extensions across hits in lanes (exact; ordering makes it hard).
4. Traceback lanes across subjects (exact; x86-64 only, likely).
5. A database index for the daemon (exact; a redesign, worth it only for many small jobs).
