package swiss.sib.blastp4j.seed;

import java.util.List;

/**
 * {@link WordFinder}'s two-hit logic and ungapped extension for words of {@link CompressedTable#WORD} letters, for
 * {@link FastBatchWordFinder} ({@code blastp -task blastp-fast}): the same code as {@link WordFinder#hit} and its
 * extensions, with the word length as a constant of its own. A copy rather than a parameter, because a word length
 * known only at run time made blastp's scan slower (c0a3894, reverted).
 * <p>
 * Not thread safe: use one per thread.
 */
final class FastWordFinder {
	private static final int WORD = CompressedTable.WORD;
	private static final byte[] PADDED = WordFinder.PADDED;
	private static final int EXTENDED = WordFinder.EXTENDED;
	private static final boolean BRANCH_FREE_MAX = WordFinder.BRANCH_FREE_MAX;

	private final byte[] query;
	private final int[] rowOfQuery;
	private final int window, xDrop, cutoff;
	private int hspQuery, hspSubject, hspLength, lastSubjectOffset, extensionLength;
	private boolean rightExtended;

	/** The finder of a query of a batch, which keeps the diagonals. */
	FastWordFinder(byte[] query, SeedParameters p) {
		this.query = query;
		this.window = p.window();
		this.xDrop = p.xDrop();
		this.cutoff = p.cutoff();
		rowOfQuery = new int[query.length];
		for (int q = 0; q < query.length; q++) {
			rowOfQuery[q] = query[q] << 5;
		}
	}

	/** As {@link WordFinder#next}. */
	int next(byte[] residues, int start, int length, int q, int s, List<UngappedHsp> out, int state, int offset) {
		if ((state & EXTENDED) != 0) {
			return s + offset >= (state & ~EXTENDED) ? s + offset : state;
		}
		final int last = state - offset;
		final int diff = s - last;
		if (diff >= window) {
			return s + offset;
		}
		if (diff < WORD) {
			return state;
		}
		if (q < diff) {
			return s + offset;
		}
		int score = extendTwoHit(residues, start, length, last + WORD, s, q);
		if (score >= cutoff) {
			out.add(new UngappedHsp(hspQuery, hspSubject, hspLength, score, q, s));
		}
		return rightExtended ? (lastSubjectOffset - (WORD - 1) + offset) | EXTENDED : s + offset;
	}

	/** As {@code WordFinder.extendTwoHit}. */
	private int extendTwoHit(byte[] residues, int start, int length, int sLeft, int sRight, int qRight) {
		int score = 0, leftScore = 0, rightD = 0;
		for (int i = 0; i < WORD; i++) {
			score += PADDED[(rowOfQuery[qRight + i] + residues[start + sRight + i]) & 1023];
			if (score > leftScore) {
				leftScore = score;
				rightD = i + 1;
			}
		}
		qRight += rightD;
		sRight += rightD;
		rightExtended = false;
		lastSubjectOffset = sRight;
		leftScore = extendLeft(residues, start, sRight - 1, qRight - 1);
		int leftD = extensionLength;
		int rightScore = 0;
		rightD = 0;
		if (leftD >= sRight - sLeft) {
			rightExtended = true;
			rightScore = extendRight(residues, start, length, sRight, qRight, leftScore);
			rightD = extensionLength;
		}
		hspQuery = qRight - leftD;
		hspSubject = sRight - leftD;
		hspLength = leftD + rightD;
		return Math.max(leftScore, rightScore);
	}

	private int extendLeft(byte[] residues, int start, int sOff, int qOff) {
		final int n = Math.min(sOff, qOff);
		int best = n + 1, score = 0, maxScore = 0, floor = -xDrop;
		final int sBase = start + sOff - n, qBase = qOff - n;
		for (int i = n; i >= 0; i--) {
			score += PADDED[(rowOfQuery[qBase + i] + residues[sBase + i]) & 1023];
			if (BRANCH_FREE_MAX) {
				best = score > maxScore ? i : best;
				maxScore = Math.max(maxScore, score);
				floor = maxScore - xDrop;
			} else if (score > maxScore) {
				maxScore = score;
				best = i;
				floor = maxScore - xDrop;
			}
			if (score <= floor) {
				break;
			}
		}
		extensionLength = n - best + 1;
		return maxScore;
	}

	private int extendRight(byte[] residues, int start, int length, int sOff, int qOff, int maxScore) {
		final int n = Math.min(length - sOff, query.length - qOff);
		int best = -1, score = maxScore, i;
		int floor = Math.max(0, maxScore - xDrop);
		for (i = 0; i < n; i++) {
			score += PADDED[(rowOfQuery[qOff + i] + residues[start + sOff + i]) & 1023];
			if (BRANCH_FREE_MAX) {
				best = score > maxScore ? i : best;
				maxScore = Math.max(maxScore, score);
				floor = Math.max(0, maxScore - xDrop);
			} else if (score > maxScore) {
				maxScore = score;
				best = i;
				floor = Math.max(0, maxScore - xDrop);
			}
			if (score <= floor) {
				break;
			}
		}
		extensionLength = best + 1;
		lastSubjectOffset = sOff + i;
		return maxScore;
	}
}
