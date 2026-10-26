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

	/** {@code ScaledSymmetricProductA}: the lower triangle of A D A^T into W. */
	private static void scaledSymmetricProductA(double[][] w, double[] diagonal, int alphsize) {
		int m = 2 * alphsize - 1;
		for (int r = 0; r < m; r++) {
			for (int c = 0; c <= r; c++) {
				w[r][c] = 0.0;
			}
		}
		for (int i = 0; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				double dd = diagonal[i * alphsize + j];
				w[j][j] += dd;
				if (i > 0) {
					w[i + alphsize - 1][j] += dd;
					w[i + alphsize - 1][i + alphsize - 1] += dd;
				}
			}
		}
	}

	/** {@code MultiplyByA}: y = beta y + alpha A x. */
	private static void multiplyByA(double beta, double[] y, int alphsize, double alpha, double[] x) {
		if (beta == 0.0) {
			for (int i = 0; i < 2 * alphsize - 1; i++) {
				y[i] = 0.0;
			}
		} else if (beta != 1.0) {
			for (int i = 0; i < 2 * alphsize - 1; i++) {
				y[i] *= beta;
			}
		}
		for (int i = 0; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				y[j] += alpha * x[i * alphsize + j];
			}
		}
		for (int i = 1; i < alphsize; i++) {
			for (int j = 0; j < alphsize; j++) {
				y[i + alphsize - 1] += alpha * x[i * alphsize + j];
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
		multiplyByA(1.0, rA, alphsize, -1.0, x);
	}

	private static void dualResiduals(double[] residsX, int alphsize, double[][] grads, double[] z,
			boolean constrainRelEntropy) {
		int n = alphsize * alphsize;
		if (constrainRelEntropy) {
			double eta = z[2 * alphsize - 1];
			for (int i = 0; i < n; i++) {
				residsX[i] = -grads[0][i] + eta * grads[1][i];
			}
		} else {
			for (int i = 0; i < n; i++) {
				residsX[i] = -grads[0][i];
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
	 * {@code Nlm_FactorLtriangPosDef}: Cholesky factorization in place of the lower triangle. Column by column, so the
	 * updates of the rows below are independent of each other. Each entry still gets its subtractions in the order of
	 * the column index, so the result is the same as the row by row factorization.
	 */
	private static void factorLtriangPosDef(double[][] a, int n) {
		final double[] column = new double[n];
		for (int k = 0; k < n; k++) {
			a[k][k] = Math.sqrt(a[k][k]);
			for (int i = k + 1; i < n; i++) {
				a[i][k] /= a[k][k];
				column[i] = a[i][k];
			}
			for (int i = k + 1; i < n; i++) {
				final double[] row = a[i];
				final double aik = row[k];
				for (int j = k + 1; j <= i; j++) {
					row[j] -= aik * column[j];
				}
			}
		}
	}

	/** {@code Nlm_SolveLtriangPosDef}: solves L L^T y = x in place. */
	private static void solveLtriangPosDef(double[] x, int n, double[][] l) {
		for (int i = 0; i < n; i++) {
			double temp = x[i];
			for (int j = 0; j < i; j++) {
				temp -= l[i][j] * x[j];
			}
			x[i] = temp / l[i][i];
		}
		for (int j = n - 1; j >= 0; j--) {
			x[j] /= l[j][j];
			for (int i = 0; i < j; i++) {
				x[i] -= l[j][i] * x[j];
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
	static int optimize(double[] x, int alphsize, double[] q, double[] rowSums, double[] colSums,
			boolean constrainRelEntropy, double relativeEntropy, double tol, int maxits) {
		final int n = alphsize * alphsize, mA = 2 * alphsize - 1, m = constrainRelEntropy ? mA + 1 : mA;
		double[] values = new double[2];
		double[][] grads = new double[2][n];
		double[] z = new double[mA + 1];
		double[] residsX = new double[n], residsZ = new double[mA + 1];
		double[] oldScores = new double[n], workspace = new double[n];
		// Newton system: W lower triangular (rows of length i + 1), Dinv, grad_re
		double[][] w = new double[2 * alphsize][];
		for (int i = 0; i < w.length; i++) {
			w[i] = new double[i + 1];
		}
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
				grads[0][k] = temp + 1;
				if (constrainRelEntropy) {
					temp += oldScores[k];
					values[1] += x[k] * temp;
					grads[1][k] = temp + 1;
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
				System.arraycopy(grads[1], 0, gradRe, 0, n);
				w[m - 1][m - 1] = 0.0;
				for (int i = 0; i < n; i++) {
					workspace[i] = dinv[i] * gradRe[i];
					w[m - 1][m - 1] += gradRe[i] * workspace[i];
				}
				multiplyByA(0.0, w[m - 1], alphsize, 1.0, workspace);
			}
			factorLtriangPosDef(w, m);
			// SolveReNewtonSystem(residsX, residsZ)
			for (int i = 0; i < n; i++) {
				workspace[i] = residsX[i] * dinv[i];
			}
			multiplyByA(1.0, residsZ, alphsize, -1.0, workspace);
			if (constrainRelEntropy) {
				for (int i = 0; i < n; i++) {
					residsZ[m - 1] -= gradRe[i] * workspace[i];
				}
			}
			solveLtriangPosDef(residsZ, m, w);
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
