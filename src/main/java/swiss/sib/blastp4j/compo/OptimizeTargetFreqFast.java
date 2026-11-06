package swiss.sib.blastp4j.compo;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import swiss.sib.blastp4j.util.WideVectors;

/**
 * {@link OptimizeTargetFreq} for {@link swiss.sib.blastp4j.util.Arithmetic#FAST} (almostblastp): the same Newton's
 * method, with every loop over the target frequencies and the Newton matrix a vector at a time. Logarithms are the
 * Vector API's (within an ulp or so of {@link Math#log}), sums are in another order, and the norms are not scaled
 * against overflow (the values are probabilities and their gradients), so the last bits differ from BLAST's.
 * Nothing reads the objective's own value, so it is not computed.
 */
final class OptimizeTargetFreqFast implements CompoKernels {
	/**
	 * The machine's preferred vectors of doubles, at most 256 bits (as the word finder's) unless
	 * {@link WideVectors#ENABLED}.
	 */
	private static final VectorSpecies<Double> D = DoubleVector.SPECIES_PREFERRED.length() > 4
			&& !WideVectors.ENABLED ? DoubleVector.SPECIES_256 : DoubleVector.SPECIES_PREFERRED;
	private static final int L = D.length();

	/** {@link OptimizeTargetFreq#freqRatioToScore}, a vector at a time. */
	@Override
	public void freqRatioToScore(double[] ratios, double lambda) {
		int k = 0;
		for (; k < D.loopBound(ratios.length); k += L) {
			DoubleVector r = DoubleVector.fromArray(D, ratios, k);
			r.lanewise(VectorOperators.LOG).div(lambda).blend(CompositionAdjustment.SCORE_MIN, r.eq(0.0))
					.intoArray(ratios, k);
		}
		for (; k < ratios.length; k++) {
			ratios[k] = 0.0 == ratios[k] ? CompositionAdjustment.SCORE_MIN : Math.log(ratios[k]) / lambda;
		}
	}

	/** {@code sum(a[ao ..] * b[bo ..])} over {@code len} entries, a vector at a time. */
	private static double dot(double[] a, int ao, double[] b, int bo, int len) {
		DoubleVector s = DoubleVector.zero(D);
		int k = 0;
		for (; k < D.loopBound(len); k += L) {
			s = DoubleVector.fromArray(D, a, ao + k).fma(DoubleVector.fromArray(D, b, bo + k), s);
		}
		double sum = s.reduceLanes(VectorOperators.ADD);
		for (; k < len; k++) {
			sum = Math.fma(a[ao + k], b[bo + k], sum);
		}
		return sum;
	}

	/** {@code y[yo ..] = beta * y[yo ..] + alpha * x[xo ..]} over {@code len} entries. */
	private static void axpby(double beta, double[] y, int yo, double alpha, double[] x, int xo, int len) {
		int k = 0;
		for (; k < D.loopBound(len); k += L) {
			DoubleVector.fromArray(D, x, xo + k).mul(alpha).add(DoubleVector.fromArray(D, y, yo + k).mul(beta))
					.intoArray(y, yo + k);
		}
		for (; k < len; k++) {
			y[yo + k] = beta * y[yo + k] + alpha * x[xo + k];
		}
	}

	/** {@code sum(v[i]^2)} over {@code n} entries; its square root is the Euclidean norm, without BLAST's scaling. */
	private static double squares(double[] v, int n) {
		return dot(v, 0, v, 0, n);
	}

	/** As {@link OptimizeTargetFreq}'s: the lower triangle of A D A^T into W (packed by rows). */
	private static void scaledSymmetricProductA(double[] w, double[] diagonal, int alphsize) {
		final int m = 2 * alphsize - 1;
		java.util.Arrays.fill(w, 0, m * (m + 1) / 2, 0.0);
		// The diagonal block: the column sums of D.
		for (int j = 0; j < alphsize; j++) {
			double s = 0;
			for (int i = 0; i < alphsize; i++) {
				s += diagonal[i * alphsize + j];
			}
			w[j * (j + 3) / 2] = s;
		}
		// Row i > 0 of D, then its sum on the diagonal, as row i + alphsize - 1.
		for (int i = 1; i < alphsize; i++) {
			final int r = (i + alphsize - 1) * (i + alphsize) / 2;
			System.arraycopy(diagonal, i * alphsize, w, r, alphsize);
			double s = 0;
			for (int j = 0; j < alphsize; j++) {
				s += diagonal[i * alphsize + j];
			}
			w[r + i + alphsize - 1] = s;
		}
	}

	/** y = beta y + alpha A x, with y at {@code y[yOff ..]}: the column sums of x, then the sums of its rows 1 on. */
	private static void multiplyByA(double beta, double[] y, int yOff, int alphsize, double alpha, double[] x) {
		final int mA = 2 * alphsize - 1;
		if (beta != 1.0) {
			for (int i = 0; i < mA; i++) {
				y[yOff + i] = beta == 0.0 ? 0.0 : beta * y[yOff + i];
			}
		}
		for (int i = 0; i < alphsize; i++) {
			axpby(1.0, y, yOff, alpha, x, i * alphsize, alphsize);
		}
		for (int i = 1; i < alphsize; i++) {
			y[yOff + i + alphsize - 1] += alpha * dot(x, i * alphsize, ONES, 0, alphsize);
		}
	}

	/** Ones, for row sums with {@link #dot}: as many as the largest alphabet. */
	private static final double[] ONES = new double[CompositionAdjustment.ALPHSIZE];
	static {
		java.util.Arrays.fill(ONES, 1.0);
	}

	/** y = beta y + alpha A^T x: entry (i, j) gets x[j], and from row 1 on x[i + alphsize - 1] too. */
	private static void multiplyByAtranspose(double beta, double[] y, int alphsize, double alpha, double[] x) {
		for (int i = 0; i < alphsize; i++) {
			final int row = i * alphsize;
			axpby(beta, y, row, alpha, x, 0, alphsize);
			if (i > 0) {
				final double add = alpha * x[i + alphsize - 1];
				int j = 0;
				for (; j < D.loopBound(alphsize); j += L) {
					DoubleVector.fromArray(D, y, row + j).add(add).intoArray(y, row + j);
				}
				for (; j < alphsize; j++) {
					y[row + j] += add;
				}
			}
		}
	}

	/**
	 * As {@link OptimizeTargetFreq}'s factorization, which leaves out the zeros of the first {@code d} rows, with each
	 * remaining dot product a vector at a time.
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
				a[ri + j] = (a[ri + j] - dot(a, ri, a, rj, j)) / a[rj + j];
			}
			a[ri + i] = Math.sqrt(a[ri + i] - dot(a, ri, a, ri, i));
		}
	}

	/** As {@link OptimizeTargetFreq}'s solve, the dot products and updates a vector at a time. */
	private static void solveLtriangPosDef(double[] x, int n, double[] l, int d) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			x[i] = (i < d ? x[i] : x[i] - dot(l, ri, x, 0, i)) / l[ri + i];
		}
		for (int j = n - 1; j >= 0; j--) {
			final int rj = j * (j + 1) / 2;
			x[j] /= l[rj + j];
			if (j >= d) {
				axpby(1.0, x, 0, -x[j], l, rj, j);
			}
		}
	}

	/** {@code Nlm_StepBound}: the largest step, at most {@code max}, that keeps x non-negative; a vector at a time. */
	private static double stepBound(double[] x, int n, double[] stepX, double max) {
		DoubleVector best = DoubleVector.broadcast(D, max);
		int i = 0;
		for (; i < D.loopBound(n); i += L) {
			DoubleVector a = DoubleVector.fromArray(D, x, i).neg().div(DoubleVector.fromArray(D, stepX, i));
			// Steps that are not >= 0 (negative, or NaN) do not bound it.
			best = best.min(a.blend(max, a.compare(VectorOperators.GE, 0.0).not()));
		}
		double alpha = best.reduceLanes(VectorOperators.MIN);
		for (; i < n; i++) {
			double alphaI = -x[i] / stepX[i];
			if (alphaI >= 0 && alphaI < alpha) {
				alpha = alphaI;
			}
		}
		return alpha;
	}

	/** {@link OptimizeTargetFreq#optimize}, with the loops a vector at a time. */
	@Override
	public int optimize(double[] x, int alphsize, double[] q, double[] rowSums, double[] colSums,
			boolean constrainRelEntropy, double relativeEntropy, double tol, int maxits) {
		final int n = alphsize * alphsize, mA = 2 * alphsize - 1, m = constrainRelEntropy ? mA + 1 : mA;
		double[] grads = new double[2 * n];
		double[] z = new double[mA + 1];
		double[] residsX = new double[n], residsZ = new double[mA + 1];
		double[] oldScores = new double[n], workspace = new double[n];
		double[] w = new double[alphsize * (2 * alphsize + 1)];
		double[] dinv = new double[n];

		for (int i = 0; i < alphsize; i++) {
			final DoubleVector r = DoubleVector.broadcast(D, rowSums[i]);
			int j = 0;
			for (; j < D.loopBound(alphsize); j += L) {
				DoubleVector.fromArray(D, q, i * alphsize + j).div(r.mul(DoubleVector.fromArray(D, colSums, j)))
						.lanewise(VectorOperators.LOG).intoArray(oldScores, i * alphsize + j);
			}
			for (; j < alphsize; j++) {
				oldScores[i * alphsize + j] = Math.log(q[i * alphsize + j] / (rowSums[i] * colSums[j]));
			}
		}
		System.arraycopy(q, 0, x, 0, n);
		int its = 0;
		double rnorm;
		while (true) {
			// EvaluateReFunctions: the gradients, and the relative entropy (the objective's value is not used).
			DoubleVector re = DoubleVector.zero(D);
			int k = 0;
			for (; k < D.loopBound(n); k += L) {
				DoubleVector xk = DoubleVector.fromArray(D, x, k);
				DoubleVector temp = xk.div(DoubleVector.fromArray(D, q, k)).lanewise(VectorOperators.LOG);
				temp.add(1).intoArray(grads, k);
				if (constrainRelEntropy) {
					temp = temp.add(DoubleVector.fromArray(D, oldScores, k));
					re = xk.fma(temp, re);
					temp.add(1).intoArray(grads, n + k);
				}
			}
			double relEntropy = re.reduceLanes(VectorOperators.ADD);
			for (; k < n; k++) {
				double temp = Math.log(x[k] / q[k]);
				grads[k] = temp + 1;
				if (constrainRelEntropy) {
					temp += oldScores[k];
					relEntropy += x[k] * temp;
					grads[n + k] = temp + 1;
				}
			}
			// CalculateResiduals
			if (constrainRelEntropy) {
				System.arraycopy(grads, n, residsX, 0, n);
				axpby(z[mA], residsX, 0, -1.0, grads, 0, n);
			} else {
				axpby(0.0, residsX, 0, -1.0, grads, 0, n);
			}
			multiplyByAtranspose(1.0, residsX, alphsize, 1.0, z);
			double normX2 = squares(residsX, n);
			for (int i = 0; i < alphsize; i++) {
				residsZ[i] = colSums[i];
			}
			for (int i = 1; i < alphsize; i++) {
				residsZ[i + alphsize - 1] = rowSums[i];
			}
			multiplyByA(1.0, residsZ, 0, alphsize, -1.0, x);
			double normZ2;
			if (constrainRelEntropy) {
				residsZ[mA] = relativeEntropy - relEntropy;
				normZ2 = squares(residsZ, mA + 1);
			} else {
				normZ2 = squares(residsZ, mA);
			}
			rnorm = Math.sqrt(normX2 + normZ2);

			if (!(rnorm > tol)) {
				break;
			}
			if (++its > maxits) {
				break;
			}
			// FactorReNewtonSystem
			axpby(0.0, dinv, 0, constrainRelEntropy ? 1 / (1 - z[m - 1]) : 1.0, x, 0, n);
			scaledSymmetricProductA(w, dinv, alphsize);
			if (constrainRelEntropy) {
				final int last = (m - 1) * m / 2;
				for (k = 0; k < D.loopBound(n); k += L) {
					DoubleVector.fromArray(D, dinv, k).mul(DoubleVector.fromArray(D, grads, n + k)).intoArray(workspace,
							k);
				}
				for (; k < n; k++) {
					workspace[k] = dinv[k] * grads[n + k];
				}
				w[last + m - 1] = dot(grads, n, workspace, 0, n);
				multiplyByA(0.0, w, last, alphsize, 1.0, workspace);
			}
			factorLtriangPosDef(w, m, alphsize);
			// SolveReNewtonSystem(residsX, residsZ)
			for (k = 0; k < D.loopBound(n); k += L) {
				DoubleVector.fromArray(D, residsX, k).mul(DoubleVector.fromArray(D, dinv, k)).intoArray(workspace, k);
			}
			for (; k < n; k++) {
				workspace[k] = residsX[k] * dinv[k];
			}
			multiplyByA(1.0, residsZ, 0, alphsize, -1.0, workspace);
			if (constrainRelEntropy) {
				residsZ[m - 1] -= dot(grads, n, workspace, 0, n);
			}
			solveLtriangPosDef(residsZ, m, w, alphsize);
			if (constrainRelEntropy) {
				axpby(1.0, residsX, 0, residsZ[m - 1], grads, n, n);
			}
			multiplyByAtranspose(1.0, residsX, alphsize, 1.0, residsZ);
			for (k = 0; k < D.loopBound(n); k += L) {
				DoubleVector.fromArray(D, residsX, k).mul(DoubleVector.fromArray(D, dinv, k)).intoArray(residsX, k);
			}
			for (; k < n; k++) {
				residsX[k] *= dinv[k];
			}
			double alpha = stepBound(x, n, residsX, 1.0 / .95) * 0.95;
			axpby(1.0, x, 0, alpha, residsX, 0, n);
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
