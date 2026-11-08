package swiss.sib.blastp4j.align;

/**
 * {@link GappedAligner#alignEx} for two extensions at a time, a row of each per {@link #step}, with the loops over
 * the two rows' cells as one: where the vectors have fewer than the 8 ints that {@link TracebackLanes} needs (AArch64's
 * 128 bits). A cell depends on the cell before it in its row (the gap along the row, then the score), a chain of
 * several cycles per cell; two rows of independent extensions give the core two chains to run side by side. Each
 * extension gives exactly alignEx's score, end point and operations: its cells are {@code ExRows.rows}'s, computed
 * with traceback throughout (without ExRows' score-only rows after the best). Not thread safe: use one per thread.
 */
public final class TracebackPairs implements TracebackKernel {
	private static final int SUB = EditScript.SUB, INS = EditScript.INS, DEL = EditScript.DEL;
	private static final int EXTEND_GAP_A = GappedAligner.EXTEND_GAP_A, EXTEND_GAP_B = GappedAligner.EXTEND_GAP_B;
	private static final int MININT = GappedAligner.MININT;

	/** One extension's alignEx state, and its row's state between the joint loop and the loop of the rest. */
	private static final class Lane {
		int[] best = new int[256], bestGap = new int[256], rowStart = new int[256];
		final int[] mat = new int[1024];
		long version = Long.MIN_VALUE;
		byte[] trace = new byte[1024], bs = new byte[256];
		byte[] a;
		int aBase, m, n, id, xDrop, extra;
		boolean reverse, busy;
		int ai, first, bSize, bestScore, aOffset, bOffset, used;
		// The current row: its matrix row, traceback base, first cell, and the cell loop's state.
		int row, base, origFirst, bi, score, gapRow, last;

		void ensureTrace(int size) {
			if (trace.length < size) {
				trace = java.util.Arrays.copyOf(trace, Math.max(size, trace.length * 2));
			}
		}
	}

	private final int goe, ext, maxA;
	private final Lane l0 = new Lane(), l1 = new Lane();
	private final EditScript ops = new EditScript();

	// The queue of submitted extensions.
	private byte[][] qa = new byte[16][], qb = new byte[16][];
	private int[] qaBase = new int[16], qm = new int[16], qbBase = new int[16], qn = new int[16], qx = new int[16],
			qid = new int[16];
	private boolean[] qrev = new boolean[16];
	private int[][] qmat = new int[16][];
	private long[] qver = new long[16];
	private int head, size;

	/** Pairs for the given gap costs, and for extensions with A of at most {@code maxA} letters. */
	public TracebackPairs(int gapOpen, int gapExtend, int maxA) {
		goe = gapOpen + gapExtend;
		ext = gapExtend;
		this.maxA = maxA;
	}

	@Override
	public boolean fits(int gapOpen, int gapExtend, int maxA) {
		return goe == gapOpen + gapExtend && ext == gapExtend && maxA <= this.maxA;
	}

	@Override
	public int maxA() {
		return maxA;
	}

	@Override
	public int queued() {
		return size;
	}

	@Override
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id,
			int[] matrix, long version) {
		if (m <= 0 || n <= 0 || m > maxA) {
			throw new IllegalArgumentException("m " + m + ", n " + n);
		}
		if (size == qa.length) {
			grow();
		}
		final int i = (head + size) % qa.length;
		qa[i] = a;
		qaBase[i] = aBase;
		qm[i] = m;
		qb[i] = b;
		qbBase[i] = bBase;
		qn[i] = n;
		qrev[i] = reverse;
		qx[i] = xDrop;
		qid[i] = id;
		qmat[i] = matrix;
		qver[i] = version;
		size++;
	}

	private void grow() {
		final int c = qa.length, c2 = 2 * c;
		final byte[][] a = new byte[c2][], b = new byte[c2][];
		final int[] aBase = new int[c2], mm = new int[c2], bBase = new int[c2], nn = new int[c2], x = new int[c2],
				ids = new int[c2];
		final boolean[] rev = new boolean[c2];
		final int[][] mat = new int[c2][];
		final long[] ver = new long[c2];
		for (int k = 0; k < size; k++) {
			final int i = (head + k) % c;
			a[k] = qa[i];
			aBase[k] = qaBase[i];
			mm[k] = qm[i];
			b[k] = qb[i];
			bBase[k] = qbBase[i];
			nn[k] = qn[i];
			rev[k] = qrev[i];
			x[k] = qx[i];
			ids[k] = qid[i];
			mat[k] = qmat[i];
			ver[k] = qver[i];
		}
		qa = a;
		qaBase = aBase;
		qm = mm;
		qb = b;
		qbBase = bBase;
		qn = nn;
		qrev = rev;
		qx = x;
		qid = ids;
		qmat = mat;
		qver = ver;
		head = 0;
	}

	@Override
	public boolean step(TracebackLanes.Done done) {
		if (!l0.busy && size > 0) {
			start(l0);
		}
		if (!l1.busy && size > 0) {
			start(l1);
		}
		if (l0.busy && l1.busy) {
			begin(l0);
			begin(l1);
			both(l0, l1);
			rest(l0);
			rest(l1);
		} else if (l0.busy || l1.busy) {
			final Lane l = l0.busy ? l0 : l1;
			begin(l);
			rest(l);
		} else {
			return false;
		}
		if (l0.busy && end(l0)) {
			finish(l0, done);
		}
		if (l1.busy && end(l1)) {
			finish(l1, done);
		}
		return true;
	}

	/** Starts the next queued extension in lane l: alignEx up to its first row. */
	private void start(Lane l) {
		final int i = head;
		head = (head + 1) % qa.length;
		size--;
		final int m = qm[i], n = qn[i], bBase = qbBase[i];
		final boolean reverse = qrev[i];
		final byte[] b = qb[i];
		l.a = qa[i];
		l.aBase = qaBase[i];
		l.m = m;
		l.n = n;
		l.reverse = reverse;
		l.id = qid[i];
		final int xDrop = Math.max(qx[i], goe);
		l.xDrop = xDrop;
		l.extra = ext > 0 ? xDrop / ext + 3 : n + 3;
		if (l.version != qver[i]) {
			l.version = qver[i];
			System.arraycopy(qmat[i], 0, l.mat, 0, 1024);
		}
		qa[i] = null;
		qb[i] = null;
		qmat[i] = null;
		if (l.best.length < n + 2) {
			final int s = Math.max(n + 2, l.best.length * 2);
			l.best = new int[s];
			l.bestGap = new int[s];
		}
		if (l.rowStart.length < m + 1) {
			l.rowStart = new int[Math.max(m + 1, l.rowStart.length * 2)];
		}
		// B's letters as GappedAligner.bLetters reads them.
		if (l.bs.length < n + 1) {
			l.bs = new byte[Math.max(n + 1, l.bs.length * 2)];
		}
		if (!reverse) {
			System.arraycopy(b, bBase + 1, l.bs, 0, Math.min(n + 1, b.length - bBase - 1));
		} else {
			for (int k = 0, end = Math.min(n + 1, bBase + n); k < end; k++) {
				l.bs[k] = b[bBase + n - 1 - k];
			}
		}
		// alignEx's first row.
		l.ensureTrace(n + 3);
		final int[] best = l.best, bestGap = l.bestGap;
		final byte[] tb = l.trace;
		l.rowStart[0] = 0;
		int score = -goe;
		best[0] = 0;
		bestGap[0] = -goe;
		int bi;
		for (bi = 1; bi <= n; bi++) {
			if (score < -xDrop) {
				break;
			}
			best[bi] = score;
			bestGap[bi] = score - goe;
			score -= ext;
			tb[bi] = (byte) DEL;
		}
		l.ai = 0;
		l.first = 0;
		l.bSize = bi;
		l.used = bi + 1;
		l.bestScore = 0;
		l.aOffset = 0;
		l.bOffset = 0;
		l.busy = true;
	}

	/** The start of lane l's next row, as ExRows.rows. */
	private void begin(Lane l) {
		final int ai = ++l.ai;
		l.ensureTrace(l.used + 2 + (l.bSize - l.first) + l.extra + 1);
		l.base = l.used + 1 - l.first;
		l.rowStart[ai] = l.base;
		l.origFirst = l.first;
		l.row = (l.reverse ? l.a[l.aBase + l.m - ai] : l.a[l.aBase + ai]) << 5;
		l.bi = l.first;
		l.score = MININT;
		l.gapRow = MININT;
		l.last = l.first;
	}

	/** The cells of both lanes' rows, as long as both have cells left: ExRows.rows' cell, once per lane. */
	private void both(Lane p, Lane q) {
		final int goe = this.goe, ext = this.ext;
		final int[] pBest = p.best, pGap = p.bestGap, pMat = p.mat, qBest = q.best, qGap = q.bestGap, qMat = q.mat;
		final byte[] pBs = p.bs, pTb = p.trace, qBs = q.bs, qTb = q.trace;
		final int pRow = p.row, pBase = p.base, pX = p.xDrop, pAi = p.ai, qRow = q.row, qBase = q.base, qX = q.xDrop,
				qAi = q.ai;
		int pBi = p.bi, pScore = p.score, pGapRow = p.gapRow, pLast = p.last, pFirst = p.first, pBestScore = p.bestScore;
		int qBi = q.bi, qScore = q.score, qGapRow = q.gapRow, qLast = q.last, qFirst = q.first, qBestScore = q.bestScore;
		final int cells = Math.min(p.bSize - pBi, q.bSize - qBi);
		for (int c = 0; c < cells; c++, pBi++, qBi++) {
			{
				int gapCol = pGap[pBi];
				final int next = pBest[pBi] + pMat[(pRow + pBs[pBi]) & 1023];
				int script = gapCol > pScore ? INS : SUB;
				final int afterCol = Math.max(pScore, gapCol);
				script = pGapRow > afterCol ? DEL : script;
				pScore = Math.max(afterCol, pGapRow);
				if (pBestScore - pScore > pX) {
					if (pFirst == pBi) {
						pFirst++;
					} else {
						pBest[pBi] = MININT;
					}
				} else {
					pLast = pBi;
					if (pScore > pBestScore) {
						pBestScore = pScore;
						p.aOffset = pAi;
						p.bOffset = pBi;
					}
					pGapRow -= ext;
					gapCol -= ext;
					if (gapCol < pScore - goe) {
						pGap[pBi] = pScore - goe;
					} else {
						pGap[pBi] = gapCol;
						script += EXTEND_GAP_B;
					}
					if (pGapRow < pScore - goe) {
						pGapRow = pScore - goe;
					} else {
						script += EXTEND_GAP_A;
					}
					pBest[pBi] = pScore;
				}
				pScore = next;
				pTb[pBase + pBi] = (byte) script;
			}
			{
				int gapCol = qGap[qBi];
				final int next = qBest[qBi] + qMat[(qRow + qBs[qBi]) & 1023];
				int script = gapCol > qScore ? INS : SUB;
				final int afterCol = Math.max(qScore, gapCol);
				script = qGapRow > afterCol ? DEL : script;
				qScore = Math.max(afterCol, qGapRow);
				if (qBestScore - qScore > qX) {
					if (qFirst == qBi) {
						qFirst++;
					} else {
						qBest[qBi] = MININT;
					}
				} else {
					qLast = qBi;
					if (qScore > qBestScore) {
						qBestScore = qScore;
						q.aOffset = qAi;
						q.bOffset = qBi;
					}
					qGapRow -= ext;
					gapCol -= ext;
					if (gapCol < qScore - goe) {
						qGap[qBi] = qScore - goe;
					} else {
						qGap[qBi] = gapCol;
						script += EXTEND_GAP_B;
					}
					if (qGapRow < qScore - goe) {
						qGapRow = qScore - goe;
					} else {
						script += EXTEND_GAP_A;
					}
					qBest[qBi] = qScore;
				}
				qScore = next;
				qTb[qBase + qBi] = (byte) script;
			}
		}
		p.bi = pBi;
		p.score = pScore;
		p.gapRow = pGapRow;
		p.last = pLast;
		p.first = pFirst;
		p.bestScore = pBestScore;
		q.bi = qBi;
		q.score = qScore;
		q.gapRow = qGapRow;
		q.last = qLast;
		q.first = qFirst;
		q.bestScore = qBestScore;
	}

	/** The cells left of lane l's row: ExRows.rows' cell. */
	private void rest(Lane l) {
		final int goe = this.goe, ext = this.ext;
		final int[] best = l.best, bestGap = l.bestGap, mat = l.mat;
		final byte[] bs = l.bs, tb = l.trace;
		final int row = l.row, base = l.base, xDrop = l.xDrop, ai = l.ai, bSize = l.bSize;
		int bi = l.bi, score = l.score, gapRow = l.gapRow, last = l.last, first = l.first, bestScore = l.bestScore;
		for (; bi < bSize; bi++) {
			int gapCol = bestGap[bi];
			final int next = best[bi] + mat[(row + bs[bi]) & 1023];
			int script = gapCol > score ? INS : SUB;
			final int afterCol = Math.max(score, gapCol);
			script = gapRow > afterCol ? DEL : script;
			score = Math.max(afterCol, gapRow);
			if (bestScore - score > xDrop) {
				if (first == bi) {
					first++;
				} else {
					best[bi] = MININT;
				}
			} else {
				last = bi;
				if (score > bestScore) {
					bestScore = score;
					l.aOffset = ai;
					l.bOffset = bi;
				}
				gapRow -= ext;
				gapCol -= ext;
				if (gapCol < score - goe) {
					bestGap[bi] = score - goe;
				} else {
					bestGap[bi] = gapCol;
					script += EXTEND_GAP_B;
				}
				if (gapRow < score - goe) {
					gapRow = score - goe;
				} else {
					script += EXTEND_GAP_A;
				}
				best[bi] = score;
			}
			score = next;
			tb[base + bi] = (byte) script;
		}
		l.bi = bi;
		l.score = score;
		l.gapRow = gapRow;
		l.last = last;
		l.first = first;
		l.bestScore = bestScore;
	}

	/** The end of lane l's row, as ExRows.rows: the next row's band. Whether the extension has ended. */
	private boolean end(Lane l) {
		if (l.first == l.bSize) {
			return true;
		}
		final int[] best = l.best, bestGap = l.bestGap;
		final int n = l.n, base = l.base;
		int bSize = l.bSize, gapRow = l.gapRow;
		if (l.last < bSize - 1) {
			bSize = l.last + 1;
		} else {
			while (gapRow >= l.bestScore - l.xDrop && bSize <= n) {
				best[bSize] = gapRow;
				bestGap[bSize] = gapRow - goe;
				gapRow -= ext;
				if (base + bSize >= l.trace.length) {
					l.ensureTrace(base + bSize + l.extra + 2);
				}
				l.trace[base + bSize] = (byte) DEL;
				bSize++;
			}
		}
		l.used += Math.max(l.bi, bSize) - l.origFirst + 1;
		if (bSize <= n) {
			best[bSize] = MININT;
			bestGap[bSize] = MININT;
			bSize++;
		}
		l.bSize = bSize;
		return l.ai >= l.m;
	}

	/** alignEx's walk from the best cell back to the start, then the result to {@code done}. */
	private void finish(Lane l, TracebackLanes.Done done) {
		ops.clear();
		final byte[] tb = l.trace;
		final int[] rowStart = l.rowStart;
		int ai = l.aOffset, bi = l.bOffset;
		int script = SUB;
		while (ai > 0 || bi > 0) {
			final int next = tb[rowStart[ai] + bi];
			switch (script) {
			case DEL:
				script = next & 0x07;
				if ((next & EXTEND_GAP_A) != 0) {
					script = DEL;
				}
				break;
			case INS:
				script = next & 0x07;
				if ((next & EXTEND_GAP_B) != 0) {
					script = INS;
				}
				break;
			default:
				script = next & 0x07;
				break;
			}
			if (script == DEL) {
				bi--;
			} else if (script == INS) {
				ai--;
			} else {
				ai--;
				bi--;
			}
			ops.add(script, 1);
		}
		l.busy = false;
		l.a = null;
		done.done(l.id, l.bestScore, l.aOffset, l.bOffset, ops);
	}
}
