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
	private final int[] diagState;
	private final int diagBase, diagMask;
	private int diagOffset;
	private int hspQuery, hspSubject, hspLength, lastSubjectOffset, extensionLength;
	private boolean rightExtended;

	/** The finder of a query at offset {@code base} of a batch whose diagonals are {@code diags}. */
	FastWordFinder(byte[] query, SeedParameters p, int[] diags, int base) {
		this.query = query;
		this.window = p.window();
		this.xDrop = p.xDrop();
		this.cutoff = p.cutoff();
		rowOfQuery = new int[query.length];
		for (int q = 0; q < query.length; q++) {
			rowOfQuery[q] = query[q] << 5;
		}
		diagState = diags;
		diagBase = base;
		diagMask = diags.length - 1;
		diagOffset = window;
	}

	/** The end of a subject of {@code length} letters: the next one's positions start beyond everything stored so far. */
	void endSubject(int length) {
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(diagState, 0);
		} else {
			diagOffset += length + window;
		}
	}

	/** As {@link WordFinder#hit}. */
	void hit(byte[] residues, int start, int length, int q, int s, List<UngappedHsp> out) {
		final int diag = (diagBase + q - s) & diagMask;
		final int state = diagState[diag];
		if ((state & EXTENDED) != 0) {
			if (s + diagOffset >= (state & ~EXTENDED)) {
				diagState[diag] = s + diagOffset;
			}
			return;
		}
		final int last = state - diagOffset;
		final int diff = s - last;
		if (diff >= window) {
			diagState[diag] = s + diagOffset;
			return;
		}
		if (diff < WORD) {
			return;
		}
		if (q < diff) {
			diagState[diag] = s + diagOffset;
			return;
		}
		int score = extendTwoHit(residues, start, length, last + WORD, s, q);
		if (score >= cutoff) {
			out.add(new UngappedHsp(hspQuery, hspSubject, hspLength, score, q, s));
		}
		if (rightExtended) {
			diagState[diag] = (lastSubjectOffset - (WORD - 1) + diagOffset) | EXTENDED;
		} else {
			diagState[diag] = s + diagOffset;
		}
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
