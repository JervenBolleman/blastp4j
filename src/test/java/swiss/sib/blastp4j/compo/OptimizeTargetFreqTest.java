package swiss.sib.blastp4j.compo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.util.Arithmetic;

/**
 * The composition stage's kernels per {@link Arithmetic} ({@link OptimizeTargetFreq#of}): BLAST's in
 * {@link OptimizeTargetFreq}, almostblastp's in {@link OptimizeTargetFreqFast}, which must reach nearly the same
 * target frequencies and scores.
 */
class OptimizeTargetFreqTest {
	static final int N = CompositionAdjustment.NUM_TRUE_AA;

	@Test
	void aClassPerArithmetic() {
		assertTrue(OptimizeTargetFreq.of(Arithmetic.BLAST) instanceof OptimizeTargetFreq);
		assertTrue(OptimizeTargetFreq.of(Arithmetic.FAST) instanceof OptimizeTargetFreqFast);
		assertSame(OptimizeTargetFreq.of(Arithmetic.FAST), OptimizeTargetFreq.of(Arithmetic.FAST));
	}

	@Test
	void theFastOptimizerConvergesAsBlasts() {
		Random random = new Random(6);
		CompoKernels blast = OptimizeTargetFreq.of(Arithmetic.BLAST), fast = OptimizeTargetFreq.of(Arithmetic.FAST);
		boolean someDiffer = false;
		for (int round = 0; round < 200; round++) {
			double[] rows = composition(random), cols = composition(random);
			boolean constrain = round % 2 == 0;
			double[] a = new double[N * N], c = new double[N * N];
			int sa = blast.optimize(a, N, Blosum62Data.JOINT_PROBS, rows, cols, constrain, 0.44, 1e-9, 2000);
			int sc = fast.optimize(c, N, Blosum62Data.JOINT_PROBS, rows, cols, constrain, 0.44, 1e-9, 2000);
			assertEquals(sa, sc, "round " + round);
			if (sa == 0) {
				// Both stop within the tolerance of the same solution.
				for (int k = 0; k < N * N; k++) {
					assertEquals(a[k], c[k], 1e-8, "round " + round + ", entry " + k);
				}
			}
			someDiffer |= !java.util.Arrays.equals(a, c);
		}
		// The fast optimizer computes another way, so some of its last bits differ.
		assertTrue(someDiffer);
	}

	@Test
	void theFastScoresAreBlasts() {
		Random random = new Random(7);
		double[] ratios = new double[CompositionAdjustment.ALPHSIZE * CompositionAdjustment.ALPHSIZE];
		for (int k = 0; k < ratios.length; k++) {
			ratios[k] = k % 7 == 0 ? 0.0 : Math.exp(random.nextGaussian());
		}
		double[] a = ratios.clone(), c = ratios.clone();
		OptimizeTargetFreq.of(Arithmetic.BLAST).freqRatioToScore(a, 0.0105);
		OptimizeTargetFreq.of(Arithmetic.FAST).freqRatioToScore(c, 0.0105);
		for (int k = 0; k < ratios.length; k++) {
			assertEquals(a[k], c[k], 1e-12 * Math.abs(a[k]), "entry " + k);
		}
	}

	/** A random composition of the 20 true amino acids, none below 1%, summing to 1. */
	private static double[] composition(Random random) {
		double[] p = new double[N];
		double sum = 0;
		for (int i = 0; i < N; i++) {
			p[i] = 0.01 + random.nextDouble();
			sum += p[i];
		}
		for (int i = 0; i < N; i++) {
			p[i] /= sum;
		}
		return p;
	}
}
