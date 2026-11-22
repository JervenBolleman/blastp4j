package swiss.sib.blastp4j.output;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

/**
 * The number formatting of {@code -outfmt 6}, at the thresholds where BLAST+ changes form. The expected strings are
 * BLAST+ 2.17.0's own output, from searches of the human reference proteome (UP000005640) against itself: no protein
 * of the bean proteome the comparison workflow uses is long enough to reach them, so
 * {@link swiss.sib.blastp4j.stats.StatisticsReferenceTest}'s 15,107 reference hits do not cover this.
 */
class TabularFormatTest {

	/**
	 * BLAST+ keeps the integer form up to 99999 and uses exponent notation above it
	 * ({@code CAlignFormatUtil::GetScoreString}). Measured: titin's self-hit (34,350 residues) prints as 70178, and
	 * self-hits of 48,000 and 50,000 residues as 97997 and 1.022e+05.
	 */
	@Test
	void bitScoreTakesExponentNotationOnlyAbove99999() {
		assertEquals("70178", TabularFormat.bitScore(70178.0));
		assertEquals("97997", TabularFormat.bitScore(97997.0));
		assertEquals("99999", TabularFormat.bitScore(99999.0));
		assertEquals("1.000e+05", TabularFormat.bitScore(100000.0));
		assertEquals("1.022e+05", TabularFormat.bitScore(102200.0));
	}

	/** A self-hit of about 5,200 residues, the shortest that reaches five figures: BLAST+ prints 10018. */
	@Test
	void bitScoreOfFiveFiguresStaysAnInteger() {
		assertEquals("10018", TabularFormat.bitScore(10018.0));
	}

	/** One decimal up to 99.9, then the integer part. */
	@Test
	void bitScoreBelowAHundredKeepsOneDecimal() {
		assertEquals("37.0", TabularFormat.bitScore(37.0));
		assertEquals("99.9", TabularFormat.bitScore(99.9));
		assertEquals("100", TabularFormat.bitScore(100.0));
		assertEquals("1177", TabularFormat.bitScore(1177.0));
	}
}
