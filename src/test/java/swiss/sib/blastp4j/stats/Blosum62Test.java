package swiss.sib.blastp4j.stats;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.cssw4j.ScoringMatrix;

class Blosum62Test {
	/** The letters of the NCBI BLOSUM62 file. */
	static final String STANDARD = "ARNDCQEGHILKMFPSTWYVBZX*";
	/** The pairs where BLAST differs from the file, as measured; X/* is -4 in both. */
	static final String DIFFERENT = "AX CX PX SX TX WX NB QZ WZ BZ";

	@Test
	void standardLettersEqualTheNcbiFile() {
		ScoringMatrix ncbi = ScoringMatrix.blosum62();
		byte[] a = ncbi.encode(STANDARD), ours = Alphabet.encode(STANDARD);
		for (int i = 0; i < STANDARD.length(); i++) {
			for (int j = 0; j < STANDARD.length(); j++) {
				char x = STANDARD.charAt(i), y = STANDARD.charAt(j);
				if (DIFFERENT.contains("" + x + y) || DIFFERENT.contains("" + y + x)) {
					continue;
				}
				assertEquals(ncbi.score(a[i], a[j]), Blosum62.score(ours[i], ours[j]),
						STANDARD.charAt(i) + "/" + STANDARD.charAt(j));
			}
		}
	}

	@Test
	void ambiguityCodePairsThatDifferFromTheFile() {
		assertEquals(4, Blosum62.score(Alphabet.encode('N'), Alphabet.encode('B')));
		assertEquals(4, Blosum62.score(Alphabet.encode('Q'), Alphabet.encode('Z')));
		assertEquals(-2, Blosum62.score(Alphabet.encode('W'), Alphabet.encode('Z')));
		assertEquals(0, Blosum62.score(Alphabet.encode('B'), Alphabet.encode('Z')));
	}

	@Test
	void xScoresMinusOneExceptAgainstStop() {
		for (int k = 1; k < Blosum62.SIZE; k++) {
			assertEquals(k == Alphabet.encode('*') ? -4 : -1, Blosum62.score(Alphabet.X, k), "X/" + Alphabet.decode(k));
		}
	}

	@Test
	void symmetric() {
		for (int i = 0; i < Blosum62.SIZE; i++) {
			for (int j = 0; j < Blosum62.SIZE; j++) {
				assertEquals(Blosum62.score(i, j), Blosum62.score(j, i), Alphabet.decode(i) + "/" + Alphabet.decode(j));
			}
		}
	}

	@Test
	void selenocysteineScoresAsCysteineAndPyrrolysineAsX() {
		byte u = Alphabet.encode('U'), c = Alphabet.encode('C'), o = Alphabet.encode('O'), x = Alphabet.X;
		for (int k = 1; k < Blosum62.SIZE; k++) {
			if (k != u && k != c && k != o && k != x) {
				assertEquals(Blosum62.score(c, k), Blosum62.score(u, k), "U/" + Alphabet.decode(k));
				assertEquals(Blosum62.score(x, k), Blosum62.score(o, k), "O/" + Alphabet.decode(k));
			}
		}
	}
}
