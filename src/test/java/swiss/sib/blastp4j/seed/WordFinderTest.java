package swiss.sib.blastp4j.seed;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.foreign.MemorySegment;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.Alphabet;

class WordFinderTest {
	static final SeedParameters PARAMS = new SeedParameters(11, 40, 16, 41, 41);

	static List<UngappedHsp> find(String query, String... subjects) {
		byte[] q = Alphabet.encode(query);
		WordFinder finder = new WordFinder(q, LookupTable.of(q, PARAMS.threshold()), PARAMS);
		List<UngappedHsp> out = new ArrayList<>();
		for (String s : subjects) {
			// In the database layout: a sentinel before and after.
			byte[] residues = new byte[s.length() + 2];
			System.arraycopy(Alphabet.encode(s), 0, residues, 1, s.length());
			finder.find(MemorySegment.ofArray(residues), 1, s.length(), out);
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
