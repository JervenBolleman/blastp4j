package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.blastp4j.stats.Blosum62;

class GappedAlignerTest {
	static final String LETTERS = "ACDEFGHIKLMNPQRSTVWY";

	/** A related pair: the subject is the query with substitutions, insertions and deletions, in database layout. */
	static byte[][] pair(Random rnd, int length) {
		StringBuilder q = new StringBuilder(), s = new StringBuilder("-");
		for (int i = 0; i < length; i++) {
			char c = LETTERS.charAt(rnd.nextInt(LETTERS.length()));
			q.append(c);
			int r = rnd.nextInt(100);
			if (r < 60) {
				s.append(c);
			} else if (r < 90) {
				s.append(LETTERS.charAt(rnd.nextInt(LETTERS.length())));
			} else if (r < 95) {
				s.append(c).append(LETTERS.charAt(rnd.nextInt(LETTERS.length())));
			}
		}
		s.append('-');
		return new byte[][] { Alphabet.encode(q), Alphabet.encode(s) };
	}

	/** The score of an edit script, recomputed from the sequences. */
	static int rescore(byte[] query, byte[] subject, int subjectBase, GappedAligner a) {
		int q = a.queryStart, s = subjectBase + a.subjectStart, score = 0;
		EditScript e = a.editScript;
		for (int i = 0; i < e.size(); i++) {
			int n = e.num(i);
			if (e.op(i) == EditScript.SUB) {
				for (int k = 0; k < n; k++) {
					score += Blosum62.score(query[q++], subject[s++]);
				}
			} else {
				score -= 11 + n;
				if (e.op(i) == EditScript.DEL) {
					s += n;
				} else {
					q += n;
				}
			}
		}
		assertEquals(a.queryStop, q, "query end");
		assertEquals(a.subjectStop, s - subjectBase, "subject end");
		return score;
	}

	@Test
	void tracebackAgreesWithScoreOnlyAndItsOwnEditScript() {
		Random rnd = new Random(42);
		GappedAligner a = new GappedAligner(11, 1);
		for (int round = 0; round < 500; round++) {
			byte[][] p = pair(rnd, 30 + rnd.nextInt(300));
			byte[] query = p[0], subject = p[1];
			int subjectLength = subject.length - 2;
			int start = rnd.nextInt(Math.min(query.length, subjectLength));
			a.extendScoreOnly(query, subject, 1, subjectLength, start, start, false, 64);
			int score = a.score, qs = a.queryStart, qe = a.queryStop, ss = a.subjectStart, se = a.subjectStop;
			a.extendWithTraceback(query, subject, 1, subjectLength, start, start, 64);
			if (a.editScript.size() > 0 && a.score == score) {
				// When both extensions end in a gap of the same kind at the start point, BLAST merges them into one
				// gap in the edit script but keeps the score of the two extensions: one gap opening more.
				int rescored = rescore(query, subject, 1, a);
				assertTrue(rescored == score || rescored == score + 11, "round " + round + ": " + rescored);
			}
			// The traceback may only drop leading or trailing gaps, which increases the score.
			assertTrue(a.score >= score, "round " + round);
			// At the last letter of either sequence BLAST's score-only extension ends the alignment one letter
			// earlier than its traceback does (found_end is false there); both are reproduced, so skip that case.
			boolean atEnd = start == query.length - 1 || start == subjectLength - 1;
			if (a.score == score && !atEnd) {
				assertEquals(qs + "-" + qe + "/" + ss + "-" + se,
						a.queryStart + "-" + a.queryStop + "/" + a.subjectStart + "-" + a.subjectStop, "round " + round);
			}
		}
	}

	@Test
	void restrictedAlignmentScoresAtMostTheExactOne() {
		Random rnd = new Random(7);
		GappedAligner a = new GappedAligner(11, 1);
		for (int round = 0; round < 500; round++) {
			byte[][] p = pair(rnd, 30 + rnd.nextInt(300));
			int subjectLength = p[1].length - 2;
			int start = rnd.nextInt(Math.min(p[0].length, subjectLength));
			a.extendScoreOnly(p[0], p[1], 1, subjectLength, start, start, false, 38);
			int exact = a.score;
			a.extendScoreOnly(p[0], p[1], 1, subjectLength, start, start, true, 38);
			assertTrue(a.score <= exact, "round " + round + ": " + a.score + " > " + exact);
		}
	}

	@Test
	void identicalSequencesAlignCompletely() {
		byte[] q = Alphabet.encode("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH");
		byte[] s = Alphabet.encode("-MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH-");
		GappedAligner a = new GappedAligner(11, 1);
		a.extendWithTraceback(q, s, 1, q.length, 20, 20, 64);
		int self = 0;
		for (byte b : q) {
			self += Blosum62.score(b, b);
		}
		assertEquals(self, a.score);
		assertEquals(q.length + "M", a.editScript.toString());
		assertEquals(0, a.queryStart);
		assertEquals(q.length, a.queryStop);
	}
}
