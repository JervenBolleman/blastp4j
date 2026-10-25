package swiss.sib.blastp4j.align;

import java.util.ArrayList;
import java.util.List;

/**
 * The HSPs already found for a subject, for the containment test that decides whether an alignment needs to be
 * computed at all: BLAST's interval tree ({@code blast_itree.c}, indexing query and subject ranges) as far as its
 * observable behaviour goes. Searching it is exact, so a plain list gives the same answers; what matters is which HSPs
 * are in it:
 * <ul>
 * <li>an HSP with the same start (query and subject offset) or the same end as one already present replaces it when
 * it is "better", and is not added otherwise ({@code s_IntervalTreeHasHSPEndpoint});</li>
 * <li>"better" is the higher score, then the shorter query range, then the shorter subject range; for identical HSPs
 * the one already present stays.</li>
 * </ul>
 * The start is checked first, so an HSP may remove one with the same start and then not be added because of one with
 * the same end, as in BLAST.
 */
public final class HspTree {
	private final List<Hsp> hsps = new ArrayList<>();

	public void clear() {
		hsps.clear();
	}

	/** {@code BlastIntervalTreeContainsHSP}: an HSP with at least this score contains both ends of the range. */
	public boolean contains(int score, int qStart, int qEnd, int sStart, int sEnd) {
		for (Hsp h : hsps) {
			if (Hsp.isContained(score, qStart, qEnd, sStart, sEnd, h)) {
				return true;
			}
		}
		return false;
	}

	/** {@code BlastIntervalTreeAddHSP} with {@code eQueryAndSubject}. */
	public void add(Hsp in) {
		if (!resolve(in, true) || !resolve(in, false)) {
			return;
		}
		hsps.add(in);
	}

	/** Resolves an HSP sharing the start (or end) point; false when the present one is better. */
	private boolean resolve(Hsp in, boolean start) {
		for (int i = 0; i < hsps.size(); i++) {
			Hsp t = hsps.get(i);
			boolean match = start ? in.queryStart == t.queryStart && in.subjectStart == t.subjectStart
					: in.queryEnd == t.queryEnd && in.subjectEnd == t.subjectEnd;
			if (!match) {
				continue;
			}
			if (keepPresent(in, t)) {
				return false;
			}
			hsps.remove(i);
			return true;
		}
		return true;
	}

	/** {@code s_HSPsHaveCommonEndpoint}: whether the HSP already present is kept rather than the new one. */
	private static boolean keepPresent(Hsp in, Hsp t) {
		if (in.score != t.score) {
			return t.score > in.score;
		}
		int inQ = in.queryEnd - in.queryStart, tQ = t.queryEnd - t.queryStart;
		if (inQ != tQ) {
			return inQ > tQ;
		}
		int inS = in.subjectEnd - in.subjectStart, tS = t.subjectEnd - t.subjectStart;
		if (inS != tS) {
			return inS > tS;
		}
		return true;
	}
}
