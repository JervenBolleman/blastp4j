package swiss.sib.blastp4j.seed;

/**
 * An ungapped alignment found by the word finder ({@code BlastInitHSP} of BLAST). Offsets are 0-based.
 *
 * @param queryStart   first query position of the alignment
 * @param subjectStart first subject position of the alignment
 * @param length       number of aligned pairs
 * @param score        raw score
 * @param querySeed    query offset of the word hit that triggered the extension
 * @param subjectSeed  subject offset of that word hit
 */
public record UngappedHsp(int queryStart, int subjectStart, int length, int score, int querySeed, int subjectSeed) {
	public int queryEnd() {
		return queryStart + length - 1;
	}

	public int subjectEnd() {
		return subjectStart + length - 1;
	}

	/** Diagonal, subject minus query offset. */
	public int diagonal() {
		return subjectStart - queryStart;
	}
}
