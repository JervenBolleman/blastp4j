package swiss.sib.blastp4j.search;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.align.Hsp;

/**
 * Merging the HSPs of two chunks of a split query ({@code Blast_HSPListsMerge} and {@code s_BlastMergeTwoHSPs}): an
 * alignment that the chunks found each in part becomes one spanning both, and the rest are kept as they are.
 */
class QueryChunkMergeTest {
	private static final int SPLIT = 1000, OVERLAP = QueryChunks.OVERLAP;

	private static Hsp hsp(int score, int qStart, int qEnd, int sStart, int sEnd) {
		return new Hsp(score, qStart, qEnd, sStart, sEnd, qStart, sStart);
	}

	@Test
	void anAlignmentAcrossTheSplitBecomesOne() {
		// The first chunk ends at 1,100, so its HSP stops there; the second starts at 1,000 and carries on. They are
		// on the same diagonal and the second one's start is inside the first, so they merge into one.
		List<Hsp> before = List.of(hsp(60, 900, 1100, 400, 600));
		List<Hsp> chunk = List.of(hsp(70, 1020, 1300, 520, 800));
		List<Hsp> merged = BlastpSearch.mergeHsps(before, chunk, SPLIT, OVERLAP);
		assertEquals(1, merged.size());
		Hsp h = merged.get(0);
		assertEquals(900, h.queryStart);
		assertEquals(1300, h.queryEnd);
		assertEquals(400, h.subjectStart);
		assertEquals(800, h.subjectEnd);
		// The score is the better of the two, or the two spread over the merged length, whichever is larger:
		// (60 + 70) / (200 + 280) * 400 = 108.
		assertEquals(108, h.score);
		assertEquals(1020, h.queryGappedStart, "the better scoring HSP's start point");
	}

	@Test
	void anAlignmentOnAnotherDiagonalIsKept() {
		// Same query range, but 40 residues off the first one's diagonal: further than OVERLAP_DIAG_CLOSE.
		List<Hsp> before = List.of(hsp(60, 900, 1100, 400, 600));
		List<Hsp> chunk = List.of(hsp(70, 1020, 1300, 560, 840));
		assertEquals(2, BlastpSearch.mergeHsps(before, chunk, SPLIT, OVERLAP).size());
	}

	@Test
	void anAlignmentThatDoesNotTouchTheOverlapIsKept() {
		// The second HSP starts past the overlap strip, so it was not cut by the split and is not a candidate.
		List<Hsp> before = List.of(hsp(60, 900, 1100, 400, 600));
		List<Hsp> chunk = List.of(hsp(70, 1200, 1400, 700, 900));
		assertEquals(2, BlastpSearch.mergeHsps(before, chunk, SPLIT, OVERLAP).size());
	}

	@Test
	void anAlignmentEndingBeforeTheSplitIsKept() {
		// The first HSP ends before the split point, so the split cannot have cut it.
		List<Hsp> before = List.of(hsp(60, 700, 900, 200, 400));
		List<Hsp> chunk = List.of(hsp(70, 1020, 1300, 520, 800));
		assertEquals(2, BlastpSearch.mergeHsps(before, chunk, SPLIT, OVERLAP).size());
	}

	@Test
	void theMergedHspsAreOrderedByScore() {
		List<Hsp> before = List.of(hsp(30, 950, 1050, 10, 110), hsp(90, 980, 1080, 500, 600));
		List<Hsp> chunk = List.of(hsp(50, 1010, 1200, 800, 990));
		List<Hsp> merged = BlastpSearch.mergeHsps(before, chunk, SPLIT, OVERLAP);
		assertEquals(3, merged.size());
		assertEquals(90, merged.get(0).score);
		assertEquals(50, merged.get(1).score);
		assertEquals(30, merged.get(2).score);
	}

	@Test
	void aChunkAtTheQueryStartIsAppended() {
		// split_offsets of 0: Blast_HitListMerge appends rather than merging, however close the diagonals are.
		List<Hsp> before = List.of(hsp(60, 900, 1100, 400, 600));
		List<Hsp> chunk = List.of(hsp(70, 1020, 1300, 520, 800));
		assertEquals(2, BlastpSearch.mergeHsps(before, chunk, 0, OVERLAP).size());
	}
}
