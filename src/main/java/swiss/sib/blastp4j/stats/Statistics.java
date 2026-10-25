package swiss.sib.blastp4j.stats;

/**
 * Bit scores and E-values as blastp reports them for gapped alignments.
 * <p>
 * The bit score is {@code (lambda S - ln K) / ln 2}. The E-value uses Spouge's finite-size correction, as BLAST+ does
 * for BLOSUM62 ({@code BLAST_SpougeStoE} in {@code blast_stat.c}): the expected number of chance alignments with score
 * at least S between the query and one subject, from the Gumbel parameters and both sequence lengths, scaled by
 * {@code databaseLength / subjectLength} to the whole database. Unlike the classic {@code m' n' K e^(-lambda S)} it
 * needs no length adjustment, and it depends on the subject length.
 * <p>
 * Instances are immutable and thread safe.
 */
public final class Statistics {
	/** 1 / sqrt(2 pi). */
	private static final double INV_SQRT_2PI = 0.39894228040143267793994605993438;

	private final ScoringParameters p;
	private final long databaseLength;
	private final double logK;
	private final double b, gumbelBeta, tau;

	/**
	 * @param parameters     the scoring system
	 * @param databaseLength total number of residues in the database
	 */
	public Statistics(ScoringParameters parameters, long databaseLength) {
		this.p = parameters;
		this.databaseLength = databaseLength;
		this.logK = Math.log(parameters.k());
		this.b = parameters.b();
		this.gumbelBeta = parameters.gumbelBeta();
		this.tau = parameters.tau();
	}

	public ScoringParameters parameters() {
		return p;
	}

	/** Bit score of a raw score. */
	public double bitScore(int rawScore) {
		return (p.lambda() * rawScore - logK) / Math.log(2);
	}

	/**
	 * E-value of a raw score for a query of {@code queryLength} residues against a subject of {@code subjectLength}
	 * residues, in a database of {@code databaseLength} residues.
	 */
	public double evalue(int rawScore, int queryLength, int subjectLength) {
		final double lambda = p.lambda();
		final double y = rawScore;
		final double a = p.alpha();
		final double alphaV = p.alphaV();
		final double sigma = p.sigma();
		// Both sides use the same (symmetric) parameters.
		final double variance = Math.max(2.0 * alphaV / lambda, alphaV * y + gumbelBeta);
		final double sqrtVariance = Math.sqrt(variance);

		final double mLength = queryLength - (a * y + b);
		final double mF = mLength / sqrtVariance;
		final double pM = Erfc.erfc(-mF / Math.sqrt(2.0)) / 2.0;
		final double p1 = mLength * pM + sqrtVariance * INV_SQRT_2PI * Math.exp(-0.5 * mF * mF);

		final double nLength = subjectLength - (a * y + b);
		final double nF = nLength / sqrtVariance;
		final double pN = Erfc.erfc(-nF / Math.sqrt(2.0)) / 2.0;
		final double p2 = nLength * pN + sqrtVariance * INV_SQRT_2PI * Math.exp(-0.5 * nF * nF);

		final double cY = Math.max(2.0 * sigma / lambda, sigma * y + tau);
		final double area = p1 * p2 + cY * pM * pN;
		return area * p.k() * Math.exp(-lambda * y) * ((double) databaseLength / subjectLength);
	}

	/**
	 * The largest raw score whose E-value is still above {@code evalue}, found by bisection as BLAST does
	 * ({@code BLAST_SpougeEtoS}); BLAST uses it as the score a hit needs to be saved.
	 */
	public int score(double evalue, int queryLength, int subjectLength) {
		int a = 0, b = Math.max((int) (Math.log(databaseLength / evalue) / p.lambda()), 2);
		double e = evalue(b, queryLength, subjectLength);
		if (e > evalue) {
			while (e > evalue) {
				a = b;
				b *= 2;
				e = evalue(b, queryLength, subjectLength);
			}
		}
		while (b - a > 1) {
			int c = (a + b) / 2;
			if (evalue(c, queryLength, subjectLength) > evalue) {
				a = c;
			} else {
				b = c;
			}
		}
		return a;
	}
}
