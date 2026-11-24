package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Random;

import org.junit.jupiter.api.Test;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * The two identities a row-vectorised X-drop kernel would stand on (step 8 of {@code docs/acceleration-plan.md}): the
 * gap along the row as a prefix maximum, and a matrix row as one table lookup. Neither is used by the code yet; these
 * tests pin the algebra so that a kernel built on it can be trusted, and they are cheap and need no data.
 */
class RowScanTest {
	/** A sentinel for "no gap here" that survives {@code + j * ext} without overflowing, unlike {@code MININT}. */
	private static final int FLOOR = -(1 << 28);

	/**
	 * {@link GappedAligner.Rows}'s row carries the gap along the row serially:
	 * {@code H(j) = max(diag(j), F(j), E(j))} with {@code E(j+1) = max(H(j) - goe, E(j) - ext)}. Because a gap that
	 * runs through an earlier {@code E} pays {@code goe} twice and {@code goe >= ext}, it is never better than opening
	 * one directly, so the maximum may be taken over the scores before {@code E} is applied:
	 * {@code E(j+1) = P(j) - goe - j * ext} with {@code P(j) = max(k <= j) (H0(k) + k * ext)} and
	 * {@code H0 = max(diag, F)}. That makes the row a prefix maximum, which is four shift-and-maximum steps over 16
	 * lanes rather than a serial chain.
	 */
	@Test
	void thePrefixMaximumFormEqualsTheSerialGapRecurrence() {
		final Random random = new Random(11);
		for (int trial = 0; trial < 4000; trial++) {
			// BLAST's costs are goe = gapOpen + gapExtend = 12 and ext = 1; the identity needs only ext <= goe.
			final int goe = 1 + random.nextInt(20), ext = 1 + random.nextInt(goe);
			final int n = 1 + random.nextInt(40);
			final int[] diag = new int[n], f = new int[n];
			for (int j = 0; j < n; j++) {
				diag[j] = random.nextInt(400) - 200;
				// Some columns have no vertical gap yet, as the band's new cells do.
				f[j] = random.nextBoolean() ? FLOOR : random.nextInt(400) - 200;
			}
			int serialE = FLOOR, prefix = FLOOR;
			for (int j = 0; j < n; j++) {
				final int h0 = Math.max(diag[j], f[j]);
				final int serial = Math.max(h0, serialE);
				serialE = Math.max(serial - goe, serialE - ext);
				final int fromPrefix = Math.max(h0, j == 0 ? FLOOR : prefix - goe - (j - 1) * ext);
				prefix = Math.max(prefix, h0 + j * ext);
				assertEquals(serial, fromPrefix, "trial " + trial + " cell " + j + " of " + n
						+ " with goe " + goe + " and ext " + ext);
			}
		}
	}

	/**
	 * A matrix row is the scores of one subject letter against all {@code BLASTAA_SIZE} letters, so 28 bytes, and the
	 * scores of BLOSUM62 fit a byte. Two table registers hold 32, so one two-vector {@code selectFrom} — {@code TBL}
	 * on AArch64, a pair of {@code vpshufb} on x86-64 — gives the scores of 16 query letters at once, in place of the
	 * gather that made the lanes lose on the Neoverse-N2 (11.07% of the program's samples in one {@code ld1w}).
	 */
	@Test
	void aMatrixRowFitsOneTableLookupOverTwoRegisters() {
		final VectorSpecies<Byte> species = ByteVector.SPECIES_128;
		final Random random = new Random(7);
		final byte[] row = new byte[32];
		for (int i = 0; i < 28; i++) {
			row[i] = (byte) (random.nextInt(16) - 4);
		}
		final ByteVector low = ByteVector.fromArray(species, row, 0);
		final ByteVector high = ByteVector.fromArray(species, row, 16);
		final byte[] letters = new byte[species.length()];
		for (int i = 0; i < letters.length; i++) {
			letters[i] = (byte) random.nextInt(28);
		}
		final ByteVector scores = ByteVector.fromArray(species, letters, 0).selectFrom(low, high);
		for (int i = 0; i < letters.length; i++) {
			assertEquals(row[letters[i]], scores.lane(i), "lane " + i + " for letter " + letters[i]);
		}
	}
}
