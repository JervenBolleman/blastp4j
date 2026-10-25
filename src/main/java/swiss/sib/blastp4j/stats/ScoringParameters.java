package swiss.sib.blastp4j.stats;

/**
 * Statistical parameters of a scoring system (substitution matrix and gap costs), as tabulated by BLAST+ in
 * {@code blast_stat.c} for the gapped case, plus the two ungapped values the finite-size correction needs.
 *
 * @param matrix    substitution matrix name
 * @param gapOpen   gap open cost (the cost of the first gap position is gapOpen + gapExtend)
 * @param gapExtend gap extension cost
 * @param lambda    Karlin-Altschul lambda
 * @param k         Karlin-Altschul K
 * @param h         relative entropy H
 * @param alpha     length-adjustment slope alpha (the Gumbel "a")
 * @param beta      length-adjustment intercept beta
 * @param c         Gumbel C
 * @param alphaV    Gumbel alpha (variance slope)
 * @param sigma     Gumbel sigma
 * @param alphaUngapped alpha of the ungapped case of the same matrix
 * @param alphaVUngapped Gumbel alpha of the ungapped case of the same matrix
 */
public record ScoringParameters(String matrix, int gapOpen, int gapExtend, double lambda, double k, double h,
		double alpha, double beta, double c, double alphaV, double sigma, double alphaUngapped, double alphaVUngapped) {

	/**
	 * BLOSUM62 with gap open 11 and extend 1, the blastp default: the {@code {11, 1, ...}} row and the ungapped row of
	 * the BLOSUM62 table in BLAST+ {@code blast_stat.c}.
	 */
	public static final ScoringParameters BLOSUM62_11_1 = new ScoringParameters("BLOSUM62", 11, 1, 0.267, 0.041, 0.14,
			1.9, -30, 0.669720, 42.602800, 43.636200, 0.7916, 4.964660);

	/** The parameters for a matrix and gap costs; only the combinations checked against BLAST+ are available. */
	public static ScoringParameters of(String matrix, int gapOpen, int gapExtend) {
		if (matrix.equalsIgnoreCase("BLOSUM62") && gapOpen == 11 && gapExtend == 1) {
			return BLOSUM62_11_1;
		}
		throw new IllegalArgumentException(
				"no statistical parameters for " + matrix + " with gap costs " + gapOpen + "/" + gapExtend);
	}

	/** Gumbel b: {@code 2 G (alphaUngapped - alpha)} with {@code G = gapOpen + gapExtend}. */
	public double b() {
		return 2.0 * (gapOpen + gapExtend) * (alphaUngapped - alpha);
	}

	/** Gumbel beta: {@code 2 G (alphaVUngapped - alphaV)}. */
	public double gumbelBeta() {
		return 2.0 * (gapOpen + gapExtend) * (alphaVUngapped - alphaV);
	}

	/** Gumbel tau: {@code 2 G (alphaVUngapped - sigma)}. */
	public double tau() {
		return 2.0 * (gapOpen + gapExtend) * (alphaVUngapped - sigma);
	}
}
