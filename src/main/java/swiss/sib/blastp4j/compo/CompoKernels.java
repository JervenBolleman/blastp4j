package swiss.sib.blastp4j.compo;

/**
 * The composition-based stage's floating-point kernels in one {@link swiss.sib.blastp4j.util.Arithmetic}: an
 * {@link OptimizeTargetFreq}, or a copy of it from {@link OptimizeTargetFreq#of} with the arithmetic as a constant.
 */
interface CompoKernels {
	/** {@link OptimizeTargetFreq#optimize}. */
	int optimize(double[] x, int alphsize, double[] q, double[] rowSums, double[] colSums, boolean constrainRelEntropy,
			double relativeEntropy, double tol, int maxits);

	/** {@link OptimizeTargetFreq#freqRatioToScore}. */
	void freqRatioToScore(double[] ratios, double lambda);
}
