package swiss.sib.blastp4j.seed;

import swiss.sib.blastp4j.stats.Statistics;
import swiss.sib.blastp4j.stats.UngappedKarlin;

/**
 * The per-query parameters of blastp's seeding stage ({@code BlastInitialWordParametersNew} and {@code ...Update} in
 * {@code blast_parameters.c}).
 *
 * @param threshold  neighbourhood word score (11)
 * @param window     two-hit window (40)
 * @param xDrop      X-drop of the ungapped extension, raw score: 7 bits with the query's ungapped lambda, rounded up
 * @param gapTrigger raw score an ungapped alignment needs for a gapped extension: 22 bits with the query's ungapped
 *                   lambda and K, rounded down
 * @param cutoff     raw score an ungapped alignment needs to be kept: the gap trigger, or the score of the E-value
 *                   cutoff if that is lower
 */
public record SeedParameters(int threshold, int window, int xDrop, int gapTrigger, int cutoff) {
	/** ln 2 as BLAST defines it ({@code NCBIMATH_LN2}). */
	static final double LN2 = 0.69314718055994530941723212145818;

	/**
	 * blastp's defaults for a query.
	 *
	 * @param query            NCBIstdaa codes
	 * @param statistics       the gapped statistics of the search
	 * @param evalue           the E-value cutoff (10); with composition-based statistics mode 2 or 3 BLAST uses 5
	 *                         times this value in the preliminary stage, so pass that
	 * @param minSubjectLength the shortest subject length for the E-value cutoff score; BLAST uses 10
	 *                         ({@link swiss.sib.blastp4j.search.BlastpSearch#MIN_SUBJECT_LENGTH})
	 */
	public static SeedParameters blastp(byte[] query, Statistics statistics, double evalue, int minSubjectLength) {
		UngappedKarlin k = UngappedKarlin.of(query);
		int xDrop = (int) Math.ceil(7 * LN2 / k.lambda());
		int gapTrigger = (int) ((22 * LN2 + k.logK()) / k.lambda());
		int cutoffMax = statistics.score(evalue, query.length, minSubjectLength);
		return new SeedParameters(11, 40, xDrop, gapTrigger, Math.min(gapTrigger, cutoffMax));
	}
}
