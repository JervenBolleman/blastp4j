package swiss.sib.blastp4j.compo;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;
import swiss.sib.blastp4j.util.Arithmetic;

/**
 * The target frequencies of a compositionally adjusted matrix ({@code optimize_target_freq.c} and
 * {@code nlm_linear_algebra.c} of BLAST+, after Yu, Wootton and Altschul 2003 and Altschul et al. 2005).
 * <p>
 * Finds the joint probabilities x closest (in Kullback-Leibler distance) to those of the standard matrix, q, whose row
 * and column sums are the compositions of the two sequences, optionally with a given relative entropy. Newton's method
 * on the Lagrangian, with the constraint matrix used implicitly, as in BLAST.
 */
final class OptimizeTargetFreq implements CompoKernels {
	private static final CompoKernels BLAST = new OptimizeTargetFreq(), FAST = new OptimizeTargetFreqFast();
	/**
	 * The machine's preferred vectors of doubles, for the loops whose outputs are each computed on their own: a vector
	 * of outputs at a time, each with BLAST's operations in BLAST's order (a multiply, then an add or subtract, never
	 * fused), so the results are BLAST's to the bit. Sums over one output stay scalar.
	 */
	private static final VectorSpecies<Double> D = DoubleVector.SPECIES_PREFERRED;
	private static final int L = D.length();

	private OptimizeTargetFreq() {
	}

	/**
	 * The composition stage's kernels in this arithmetic: BLAST's (this class), or almostblastp's
	 * ({@link OptimizeTargetFreqFast}, a copy with its loops a vector at a time).
	 */
	static CompoKernels of(Arithmetic arithmetic) {
		return arithmetic == Arithmetic.FAST ? FAST : BLAST;
	}

	/** {@code Blast_FreqRatioToScore} on every entry of {@code ratios}. */
	@Override
	public void freqRatioToScore(double[] ratios, double lambda) {
		for (int k = 0; k < ratios.length; k++) {
			ratios[k] = 0.0 == ratios[k] ? CompositionAdjustment.SCORE_MIN : Math.log(ratios[k]) / lambda;
		}
	}

	// A lower triangle packed by rows: row i starts at i * (i + 1) / 2. Written out where used, not as a method:
	// bin/jvm.options keeps C2 from inlining this class's methods.

	/** {@code ScaledSymmetricProductA}: the lower triangle of A D A^T into W (packed by rows). */
	private static void scaledSymmetricProductA(double[] w, double[] diagonal, int alphsize) {
		int m = 2 * alphsize - 1;
		java.util.Arrays.fill(w, 0, m * (m + 1) / 2, 0.0);
		for (int i = 0; i < alphsize; i++) {
			final int r = (i + alphsize - 1) * (i + alphsize) / 2;
			for (int j = 0; j < alphsize; j++) {
				double dd = diagonal[i * alphsize + j];
				w[j * (j + 3) / 2] += dd;
				if (i > 0) {
					w[r + j] += dd;
					w[r + i + alphsize - 1] += dd;
				}
			}
		}
	}

	/** {@code MultiplyByA}: y = beta y + alpha A x, with y at {@code y[yOff ..]}. */
	private static void multiplyByA(double beta, double[] y, int yOff, int alphsize, double alpha, double[] x) {
		if (beta == 0.0) {
			for (int i = 0; i < 2 * alphsize - 1; i++) {
				y[yOff + i] = 0.0;
			}
		} else if (beta != 1.0) {
			for (int i = 0; i < 2 * alphsize - 1; i++) {
				y[yOff + i] *= beta;
			}
		}
		final int bound = D.loopBound(alphsize);
		for (int i = 0; i < alphsize; i++) {
			int j = 0;
			for (; j < bound; j += L) {
				DoubleVector.fromArray(D, y, yOff + j)
						.add(DoubleVector.fromArray(D, x, i * alphsize + j).mul(alpha)).intoArray(y, yOff + j);
			}
			for (; j < alphsize; j++) {
				y[yOff + j] += alpha * x[i * alphsize + j];
			}
		}
		for (int i = 1; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				y[yOff + i + alphsize - 1] += alpha * x[i * alphsize + j];
			}
		}
	}

	/** {@code MultiplyByAtranspose}: y = beta y + alpha A^T x. */
	private static void multiplyByAtranspose(double beta, double[] y, int alphsize, double alpha, double[] x) {
		int n = alphsize * alphsize;
		if (beta == 0.0) {
			for (int k = 0; k < n; k++) {
				y[k] = 0.0;
			}
		} else if (beta != 1.0) {
			for (int k = 0; k < n; k++) {
				y[k] *= beta;
			}
		}
		final int bound = D.loopBound(alphsize);
		for (int i = 0; i < alphsize; i++) {
			int j = 0;
			if (i > 0) {
				final double c = alpha * x[i + alphsize - 1];
				for (; j < bound; j += L) {
					final int k = i * alphsize + j;
					DoubleVector.fromArray(D, y, k).add(DoubleVector.fromArray(D, x, j).mul(alpha)).add(c)
							.intoArray(y, k);
				}
				for (; j < alphsize; j++) {
					final int k = i * alphsize + j;
					y[k] += alpha * x[j];
					y[k] += c;
				}
			} else {
				for (; j < bound; j += L) {
					DoubleVector.fromArray(D, y, j).add(DoubleVector.fromArray(D, x, j).mul(alpha)).intoArray(y, j);
				}
				for (; j < alphsize; j++) {
					y[j] += alpha * x[j];
				}
			}
		}
	}

	private static void residualsLinearConstraints(double[] rA, int alphsize, double[] x, double[] rowSums,
			double[] colSums) {
		for (int i = 0; i < alphsize; i++) {
			rA[i] = colSums[i];
		}
		for (int i = 1; i < alphsize; i++) {
			rA[i + alphsize - 1] = rowSums[i];
		}
		multiplyByA(1.0, rA, 0, alphsize, -1.0, x);
	}

	/** {@code grads}: the two gradients, of n entries each, one after the other. */
	private static void dualResiduals(double[] residsX, int alphsize, double[] grads, double[] z,
			boolean constrainRelEntropy) {
		int n = alphsize * alphsize;
		if (constrainRelEntropy) {
			double eta = z[2 * alphsize - 1];
			int i = 0;
			for (; i < D.loopBound(n); i += L) {
				DoubleVector.fromArray(D, grads, i).neg().add(DoubleVector.fromArray(D, grads, n + i).mul(eta))
						.intoArray(residsX, i);
			}
			for (; i < n; i++) {
				residsX[i] = -grads[i] + eta * grads[n + i];
			}
		} else {
			for (int i = 0; i < n; i++) {
				residsX[i] = -grads[i];
			}
		}
		multiplyByAtranspose(1.0, residsX, alphsize, 1.0, z);
	}

	/** {@code Nlm_EuclideanNorm}, scaled to avoid overflow. */
	static double euclideanNorm(double[] v, int n) {
		double sum = 1.0, scale = 0.0;
		for (int i = 0; i < n; i++) {
			if (v[i] != 0.0) {
				double absvi = Math.abs(v[i]);
				if (scale < absvi) {
					sum = 1.0 + sum * (scale / absvi) * (scale / absvi);
					scale = absvi;
				} else {
					sum += (absvi / scale) * (absvi / scale);
				}
			}
		}
		return scale * Math.sqrt(sum);
	}

	/**
	 * {@code Nlm_FactorLtriangPosDef}: Cholesky factorization in place of the lower triangle (packed by rows), whose
	 * first {@code d} rows are diagonal (as {@link #scaledSymmetricProductA} makes them: A D A^T's block of the column
	 * constraints). BLAST's loops run over those zeros too; here they are left out. Products with them are zeros,
	 * whose subtraction leaves every other value as it is, and the remaining terms are summed in BLAST's order, so the
	 * factor is BLAST's: in the first {@code d} rows only the square roots of the diagonal, and in the others the first
	 * {@code d} columns are only divided by those.
	 */
	private static void factorLtriangPosDef(double[] a, int n, int d) {
		for (int i = 0; i < Math.min(d, n); i++) {
			final int ri = i * (i + 1) / 2;
			a[ri + i] = Math.sqrt(a[ri + i]);
		}
		for (int i = d; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			for (int j = 0; j < d; j++) {
				a[ri + j] = a[ri + j] / a[j * (j + 1) / 2 + j];
			}
			for (int j = d; j < i; j++) {
				final int rj = j * (j + 1) / 2;
				double temp = a[ri + j];
				for (int k = 0; k < j; k++) {
					temp -= a[ri + k] * a[rj + k];
				}
				a[ri + j] = temp / a[rj + j];
			}
			double temp = a[ri + i];
			for (int k = 0; k < i; k++) {
				temp -= a[ri + k] * a[ri + k];
			}
			a[ri + i] = Math.sqrt(temp);
		}
	}

	/**
	 * {@code Nlm_SolveLtriangPosDef}: solves L L^T y = x in place; L packed by rows, its first {@code d} rows diagonal
	 * (as {@link #factorLtriangPosDef} leaves them), whose zeros are left out as there.
	 */
	private static void solveLtriangPosDef(double[] x, int n, double[] l, int d) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			double temp = x[i];
			for (int j = 0; j < (i < d ? 0 : i); j++) {
				temp -= l[ri + j] * x[j];
			}
			x[i] = temp / l[ri + i];
		}
		for (int j = n - 1; j >= 0; j--) {
			final int rj = j * (j + 1) / 2;
			x[j] /= l[rj + j];
			final int end = j < d ? 0 : j;
			final double xj = x[j];
			int i = 0;
			for (; i < D.loopBound(end); i += L) {
				DoubleVector.fromArray(D, x, i).sub(DoubleVector.fromArray(D, l, rj + i).mul(xj)).intoArray(x, i);
			}
			for (; i < end; i++) {
				x[i] -= l[rj + i] * xj;
			}
		}
	}

	/** {@code Nlm_StepBound}. */
	private static double stepBound(double[] x, int n, double[] stepX, double max) {
		double alpha = max;
		for (int i = 0; i < n; i++) {
			double alphaI = -x[i] / stepX[i];
			if (alphaI >= 0 && alphaI < alpha) {
				alpha = alphaI;
			}
		}
		return alpha;
	}

	/**
	 * {@code Blast_OptimizeTargetFrequencies}: fills {@code x} (row major, alphsize x alphsize); returns 0 when Newton's
	 * method converged to a minimizer, 1 otherwise.
	 */
	@Override
	public int optimize(double[] x, int alphsize, double[] q, double[] rowSums, double[] colSums,
			boolean constrainRelEntropy, double relativeEntropy, double tol, int maxits) {
		final int n = alphsize * alphsize, mA = 2 * alphsize - 1, m = constrainRelEntropy ? mA + 1 : mA;
		double[] values = new double[2];
		// The gradients of the two functions, of n entries each, one after the other.
		double[] grads = new double[2 * n];
		double[] z = new double[mA + 1];
		double[] residsX = new double[n], residsZ = new double[mA + 1];
		double[] oldScores = new double[n], workspace = new double[n];
		// Newton system: W lower triangular, packed by rows (row i has i + 1 entries), Dinv, grad_re
		double[] w = new double[alphsize * (2 * alphsize + 1)];
		double[] dinv = new double[n], gradRe = new double[n];

		for (int i = 0; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				int k = i * alphsize + j;
				oldScores[k] = Math.log(q[k] / (rowSums[i] * colSums[j]));
			}
		}
		System.arraycopy(q, 0, x, 0, n);
		int its = 0;
		double rnorm;
		while (true) {
			// EvaluateReFunctions
			values[0] = 0.0;
			values[1] = 0.0;
			for (int k = 0; k < n; k++) {
				double temp = Math.log(x[k] / q[k]);
				values[0] += x[k] * temp;
				grads[k] = temp + 1;
				if (constrainRelEntropy) {
					temp += oldScores[k];
					values[1] += x[k] * temp;
					grads[n + k] = temp + 1;
				}
			}
			// CalculateResiduals
			dualResiduals(residsX, alphsize, grads, z, constrainRelEntropy);
			double normX = euclideanNorm(residsX, n);
			residualsLinearConstraints(residsZ, alphsize, x, rowSums, colSums);
			double normZ;
			if (constrainRelEntropy) {
				residsZ[2 * alphsize - 1] = relativeEntropy - values[1];
				normZ = euclideanNorm(residsZ, 2 * alphsize);
			} else {
				normZ = euclideanNorm(residsZ, 2 * alphsize - 1);
			}
			rnorm = Math.sqrt(normX * normX + normZ * normZ);

			if (!(rnorm > tol)) {
				break;
			}
			if (++its > maxits) {
				break;
			}
			// FactorReNewtonSystem
			if (constrainRelEntropy) {
				double eta = z[m - 1];
				for (int i = 0; i < n; i++) {
					dinv[i] = x[i] / (1 - eta);
				}
			} else {
				System.arraycopy(x, 0, dinv, 0, n);
			}
			scaledSymmetricProductA(w, dinv, alphsize);
			if (constrainRelEntropy) {
				System.arraycopy(grads, n, gradRe, 0, n);
				final int last = (m - 1) * m / 2;
				w[last + m - 1] = 0.0;
				for (int i = 0; i < n; i++) {
					workspace[i] = dinv[i] * gradRe[i];
					w[last + m - 1] += gradRe[i] * workspace[i];
				}
				multiplyByA(0.0, w, last, alphsize, 1.0, workspace);
			}
			// The first alphsize rows of W are diagonal (the column constraints).
			factorLtriangPosDef(w, m, alphsize);
			// SolveReNewtonSystem(residsX, residsZ)
			for (int i = 0; i < n; i++) {
				workspace[i] = residsX[i] * dinv[i];
			}
			multiplyByA(1.0, residsZ, 0, alphsize, -1.0, workspace);
			if (constrainRelEntropy) {
				for (int i = 0; i < n; i++) {
					residsZ[m - 1] -= gradRe[i] * workspace[i];
				}
			}
			solveLtriangPosDef(residsZ, m, w, alphsize);
			if (constrainRelEntropy) {
				for (int i = 0; i < n; i++) {
					residsX[i] += gradRe[i] * residsZ[m - 1];
				}
			}
			multiplyByAtranspose(1.0, residsX, alphsize, 1.0, residsZ);
			for (int i = 0; i < n; i++) {
				residsX[i] *= dinv[i];
			}
			double alpha = stepBound(x, n, residsX, 1.0 / .95) * 0.95;
			for (int i = 0; i < n; i++) {
				x[i] += alpha * residsX[i];
			}
			for (int i = 0; i < m; i++) {
				z[i] += alpha * residsZ[i];
			}
		}
		boolean converged = false;
		if (its <= maxits && rnorm <= tol) {
			if (!constrainRelEntropy || z[m - 1] < 1) {
				converged = true;
			}
		}
		return converged ? 0 : 1;
	}
}
