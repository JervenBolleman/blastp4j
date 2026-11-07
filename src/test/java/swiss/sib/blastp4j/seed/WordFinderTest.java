package swiss.sib.blastp4j.seed;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.Alphabet;

class WordFinderTest {
	static final SeedParameters PARAMS = new SeedParameters(11, 40, 16, 41, 41);

	/**
	 * The alignments of both scans, {@link WordFinder#find} and {@link WordFinder#findBlocks}, of a word finder and of
	 * its copy with constants ({@link WordFinder#of}), which must all agree.
	 */
	static List<UngappedHsp> find(String query, String... subjects) {
		List<UngappedHsp> one = find(false, false, query, subjects);
		assertEquals(one, find(true, false, query, subjects));
		assertEquals(one, find(false, true, query, subjects));
		assertEquals(one, find(true, true, query, subjects));
		return one;
	}

	static List<UngappedHsp> find(boolean blocks, boolean copy, String query, String... subjects) {
		byte[] q = Alphabet.encode(query);
		LookupTable table = LookupTable.of(q, PARAMS.threshold());
		WordScan finder = copy ? WordFinder.of(q, table, PARAMS) : new WordFinder(q, table, PARAMS);
		List<UngappedHsp> out = new ArrayList<>();
		for (String s : subjects) {
			// In the database layout: a sentinel before and after.
			byte[] residues = new byte[s.length() + 2];
			System.arraycopy(Alphabet.encode(s), 0, residues, 1, s.length());
			if (blocks) {
				finder.findBlocks(residues, 1, s.length(), out);
			} else {
				finder.find(residues, 1, s.length(), out);
			}
		}
		return out;
	}

	static int selfScore(String s) {
		byte[] e = Alphabet.encode(s);
		int score = 0;
		for (byte b : e) {
			score += swiss.sib.blastp4j.stats.Blosum62.score(b, b);
		}
		return score;
	}

	static final String HBA = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";

	@Test
	void aCopyPerSetOfConstants() {
		byte[] q = Alphabet.encode(HBA);
		LookupTable table = LookupTable.of(q, PARAMS.threshold());
		WordScan a = WordFinder.of(q, table, PARAMS), b = WordFinder.of(q, table, PARAMS);
		assertTrue(a.getClass().isHidden(), a.getClass().toString());
		assertEquals(a.getClass(), b.getClass());
		SeedParameters other = new SeedParameters(11, 40, 15, 41, 41);
		assertNotEquals(a.getClass(), WordFinder.of(q, table, other).getClass());
	}

	@Test
	void identicalSequencesGiveOneFullLengthAlignment() {
		List<UngappedHsp> hsps = find(HBA, HBA);
		assertEquals(1, hsps.size(), hsps.toString());
		UngappedHsp h = hsps.get(0);
		assertEquals(0, h.queryStart());
		assertEquals(0, h.subjectStart());
		assertEquals(HBA.length(), h.length());
		assertEquals(selfScore(HBA), h.score());
	}

	@Test
	void shiftedCopyIsOnItsDiagonal() {
		String subject = "GGGGGGGGGG" + HBA.substring(20, 80) + "PPPPPPPPPP";
		List<UngappedHsp> hsps = find(HBA, subject);
		assertEquals(1, hsps.size(), hsps.toString());
		UngappedHsp h = hsps.get(0);
		assertEquals(10 - 20, h.diagonal());
		assertTrue(h.queryStart() <= 20 && h.queryEnd() >= 79, h.toString());
	}

	@Test
	void aSingleWordHitIsNotExtended() {
		// One shared word (WGK), nothing else in common: no second hit on the diagonal.
		assertTrue(find(HBA, "PPPPPPPPPPWGKPPPPPPPPPP").isEmpty());
	}

	@Test
	void stateDoesNotCarryOverToTheNextSubject() {
		// The same subject twice gives the same alignment twice.
		List<UngappedHsp> hsps = find(HBA, HBA, HBA);
		assertEquals(2, hsps.size());
		assertEquals(hsps.get(0), hsps.get(1));
	}
}
