package swiss.sib.blastp4j.compo;

/**
 * The target frequencies of a compositionally adjusted matrix ({@code optimize_target_freq.c} and
 * {@code nlm_linear_algebra.c} of BLAST+, after Yu, Wootton and Altschul 2003 and Altschul et al. 2005).
 * <p>
 * Finds the joint probabilities x closest (in Kullback-Leibler distance) to those of the standard matrix, q, whose row
 * and column sums are the compositions of the two sequences, optionally with a given relative entropy. Newton's method
 * on the Lagrangian, with the constraint matrix used implicitly, as in BLAST.
 */
final class OptimizeTargetFreq {
	private OptimizeTargetFreq() {
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
		for (int i = 0; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
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
		for (int i = 0; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				int k = i * alphsize + j;
				y[k] += alpha * x[j];
				if (i > 0) {
					y[k] += alpha * x[i + alphsize - 1];
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
			for (int i = 0; i < n; i++) {
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

	/** {@code Nlm_FactorLtriangPosDef}: Cholesky factorization in place of the lower triangle (packed by rows). */
	private static void factorLtriangPosDef(double[] a, int n) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			for (int j = 0; j < i; j++) {
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
	 * {@link #factorLtriangPosDef} for {@link swiss.sib.blastp4j.util.Arithmetic#FAST}: each dot product in four
	 * partial sums with fused multiply-adds, so that they do not wait for each other, and a multiplication by the
	 * diagonal entry's reciprocal ({@code inverse}, n entries, filled here) instead of a division: each entry of a row
	 * waits for the one before it, through the division. The sums are in another order and the quotients rounded
	 * twice, so the last bits can differ from BLAST's.
	 */
	private static void factorLtriangPosDefFast(double[] a, int n, double[] inverse) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			for (int j = 0; j < i; j++) {
				final int rj = j * (j + 1) / 2;
				double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
				int k = 0;
				for (; k + 4 <= j; k += 4) {
					s0 = Math.fma(a[ri + k], a[rj + k], s0);
					s1 = Math.fma(a[ri + k + 1], a[rj + k + 1], s1);
					s2 = Math.fma(a[ri + k + 2], a[rj + k + 2], s2);
					s3 = Math.fma(a[ri + k + 3], a[rj + k + 3], s3);
				}
				for (; k < j; k++) {
					s0 = Math.fma(a[ri + k], a[rj + k], s0);
				}
				a[ri + j] = (a[ri + j] - ((s0 + s1) + (s2 + s3))) * inverse[j];
			}
			double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
			int k = 0;
			for (; k + 4 <= i; k += 4) {
				s0 = Math.fma(a[ri + k], a[ri + k], s0);
				s1 = Math.fma(a[ri + k + 1], a[ri + k + 1], s1);
				s2 = Math.fma(a[ri + k + 2], a[ri + k + 2], s2);
				s3 = Math.fma(a[ri + k + 3], a[ri + k + 3], s3);
			}
			for (; k < i; k++) {
				s0 = Math.fma(a[ri + k], a[ri + k], s0);
			}
			a[ri + i] = Math.sqrt(a[ri + i] - ((s0 + s1) + (s2 + s3)));
			inverse[i] = 1 / a[ri + i];
		}
	}

	/** {@code Nlm_SolveLtriangPosDef}: solves L L^T y = x in place; L packed by rows. */
	private static void solveLtriangPosDef(double[] x, int n, double[] l) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			double temp = x[i];
			for (int j = 0; j < i; j++) {
				temp -= l[ri + j] * x[j];
			}
			x[i] = temp / l[ri + i];
		}
		for (int j = n - 1; j >= 0; j--) {
			final int rj = j * (j + 1) / 2;
			x[j] /= l[rj + j];
			for (int i = 0; i < j; i++) {
				x[i] -= l[rj + i] * x[j];
			}
		}
	}

	/**
	 * {@link #solveLtriangPosDef} for {@link swiss.sib.blastp4j.util.Arithmetic#FAST}, as
	 * {@link #factorLtriangPosDefFast}.
	 */
	private static void solveLtriangPosDefFast(double[] x, int n, double[] l) {
		for (int i = 0; i < n; i++) {
			final int ri = i * (i + 1) / 2;
			double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
			int j = 0;
			for (; j + 4 <= i; j += 4) {
				s0 = Math.fma(l[ri + j], x[j], s0);
				s1 = Math.fma(l[ri + j + 1], x[j + 1], s1);
				s2 = Math.fma(l[ri + j + 2], x[j + 2], s2);
				s3 = Math.fma(l[ri + j + 3], x[j + 3], s3);
			}
			for (; j < i; j++) {
				s0 = Math.fma(l[ri + j], x[j], s0);
			}
			x[i] = (x[i] - ((s0 + s1) + (s2 + s3))) / l[ri + i];
		}
		for (int j = n - 1; j >= 0; j--) {
			final int rj = j * (j + 1) / 2;
			x[j] /= l[rj + j];
			final double xj = x[j];
			for (int i = 0; i < j; i++) {
				x[i] = Math.fma(-l[rj + i], xj, x[i]);
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
	 * method converged to a minimizer, 1 otherwise. With {@code fast}, the Newton system is factored and solved
	 * with {@link #factorLtriangPosDefFast} and {@link #solveLtriangPosDefFast}.
	 */
	static int optimize(double[] x, int alphsize, double[] q, double[] rowSums, double[] colSums,
			boolean constrainRelEntropy, double relativeEntropy, double tol, int maxits, boolean fast) {
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
		// The reciprocals of the factor's diagonal, for the fast factorization.
		double[] inverse = new double[m];

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
			if (fast) {
				factorLtriangPosDefFast(w, m, inverse);
			} else {
				factorLtriangPosDef(w, m);
			}
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
			if (fast) {
				solveLtriangPosDefFast(residsZ, m, w);
			} else {
				solveLtriangPosDef(residsZ, m, w);
			}
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
