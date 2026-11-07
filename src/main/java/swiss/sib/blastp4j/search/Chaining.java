package swiss.sib.blastp4j.search;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import swiss.sib.blastp4j.seed.UngappedHsp;

/**
 * {@code s_ChainingAlignment} ({@code blast_gapalign.c}), which {@code blastp -task blastp-fast} runs before the
 * preliminary gapped extensions: it estimates the gapped score each ungapped alignment could reach by chaining it with
 * co-linear ungapped alignments after it, and drops those whose estimate, plus a gap and the ungapped cutoff, stays
 * below the gapped cutoff.
 * <p>
 * BLAST runs it on a subject's ungapped alignments of all queries in the batch, grouped by query; here it runs on one
 * query's. BLAST then removes the dropped ones by moving the last alignment of the batch into each gap, and sorts by
 * score; that order only differs from sorting the query's alignments where two have the same score, start and length
 * (none so far).
 */
final class Chaining {
	/** {@code s_CompareInitHSPsByQueryOffsetScore}: by query start, then best score first. */
	private static final Comparator<UngappedHsp> BY_QUERY_START = (h1, h2) -> {
		int r = Integer.compare(h1.queryStart(), h2.queryStart());
		return r != 0 ? r : Integer.compare(h2.score(), h1.score());
	};

	private Chaining() {
	}

	/**
	 * Drops from {@code seeds} (one query's ungapped alignments in one subject) those the chaining estimate rules out,
	 * leaving the others in BLAST's order by score ({@link BlastpSearch#SEED_ORDER}).
	 *
	 * @param ungappedCutoff the query's ungapped cutoff ({@code word_params->cutoffs[context].cutoff_score})
	 * @param gappedCutoff   the query's gapped cutoff ({@code hit_params->cutoffs[context].cutoff_score})
	 */
	static void filter(List<UngappedHsp> seeds, int ungappedCutoff, int gappedCutoff, int gapOpen, int gapExtend) {
		// The word finder's output is sorted by score first (BlastAaWordFinder), then by query start (stable sorts).
		seeds.sort(BlastpSearch.SEED_ORDER);
		seeds.sort(BY_QUERY_START);
		final int n = seeds.size();
		final int[] best = new int[n];
		for (int k = 0; k < n; k++) {
			best[k] = seeds.get(k).score();
		}
		for (int k = n - 1; k >= 0; k--) {
			final UngappedHsp hk = seeds.get(k);
			final int self = best[k];
			for (int j = k + 1; j < n; j++) {
				final UngappedHsp hj = seeds.get(j);
				// As BLAST: the distance from k's start, plus k's length.
				final int qDiff = hj.queryStart() - hk.queryStart() + hk.length();
				final int sDiff = hj.subjectStart() - hk.subjectStart() + hk.length();
				if (sDiff < 0) {
					continue; // only co-linear alignments
				}
				int score = self + best[j];
				// Between the two, a few matches scoring less than the ungapped cutoff, and at least one gap.
				score += Math.min(Math.min(qDiff, sDiff) * 3, ungappedCutoff);
				score -= Math.max(Math.abs(qDiff - sDiff), 1) + gapOpen;
				if (score > best[k]) {
					best[k] = score;
				}
			}
		}
		final int gapScore = gapOpen + gapExtend;
		final List<UngappedHsp> kept = new ArrayList<>(n);
		for (int k = 0; k < n; k++) {
			if (best[k] - gapScore + ungappedCutoff - 1 >= gappedCutoff) {
				kept.add(seeds.get(k));
			}
		}
		kept.sort(BlastpSearch.SEED_ORDER);
		seeds.clear();
		seeds.addAll(kept);
	}
}
