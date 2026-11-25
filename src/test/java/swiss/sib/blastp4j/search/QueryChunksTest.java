package swiss.sib.blastp4j.search;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * BLAST's chunks of a long query batch. The thresholds are measured against BLAST+ 2.17.0 rather than read from its
 * source: with titin (34,350 residues) against one subject, the HSP at query 18,869-18,915 appears exactly when the
 * query reaches 19,800 residues, which is {@code 2 * (10000 - 100)}, and with {@code OVERLAP_CHUNK_SIZE} 200 and 500
 * it appears at 19,600 and 19,000, which are {@code 2 * (10000 - 200)} and {@code 2 * (10000 - 500)}. So the chunk
 * size is 10,000 and the number of chunks is the integer division by the chunk size less the overlap.
 */
class QueryChunksTest {
	@Test
	void aQueryShorterThanTwoChunksIsNotSplit() {
		for (int length : new int[] { 1, 10000, 19799 }) {
			QueryChunks c = QueryChunks.of(length, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP);
			assertFalse(c.split(), length + " residues");
			assertEquals(1, c.count());
			assertEquals(0, c.start(0));
			assertEquals(length, c.end(0));
		}
	}

	@Test
	void theSplitStartsAtTwiceTheChunkSizeLessTheOverlap() {
		assertFalse(QueryChunks.of(19799, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP).split());
		assertTrue(QueryChunks.of(19800, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP).split());
		assertFalse(QueryChunks.of(19599, QueryChunks.CHUNK_SIZE, 200).split());
		assertTrue(QueryChunks.of(19600, QueryChunks.CHUNK_SIZE, 200).split());
		assertFalse(QueryChunks.of(18999, QueryChunks.CHUNK_SIZE, 500).split());
		assertTrue(QueryChunks.of(19000, QueryChunks.CHUNK_SIZE, 500).split());
	}

	@Test
	void twoChunksEvenedOutAndOverlapping() {
		// 19,800 / 9,900 = 2 chunks; the size is (19800 + 100) / 2 = 9,950, rounded up to 9,951 since 2 < 9,851.
		QueryChunks c = QueryChunks.of(19800, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP);
		assertEquals(2, c.count());
		assertEquals(0, c.start(0));
		assertEquals(9951, c.end(0));
		assertEquals(9851, c.start(1));
		assertEquals(19800, c.end(1));
		assertEquals(QueryChunks.OVERLAP, c.end(0) - c.start(1));
	}

	@Test
	void titinIsThreeChunks() {
		// 34,350 / 9,900 = 3 chunks; the size is (34350 + 200) / 3 = 11,516, rounded up to 11,517.
		QueryChunks c = QueryChunks.of(34350, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP);
		assertEquals(3, c.count());
		assertEquals(0, c.start(0));
		assertEquals(11517, c.end(0));
		assertEquals(11417, c.start(1));
		assertEquals(22934, c.end(1));
		assertEquals(22834, c.start(2));
		assertEquals(34350, c.end(2));
	}

	@Test
	void theLastChunkReachesTheEnd() {
		for (int length = 19800; length < 60000; length += 97) {
			QueryChunks c = QueryChunks.of(length, QueryChunks.CHUNK_SIZE, QueryChunks.OVERLAP);
			assertEquals(length, c.end(c.count() - 1), length + " residues");
			assertEquals(0, c.start(0));
			for (int i = 1; i < c.count(); i++) {
				assertTrue(c.start(i) < c.end(i - 1), "chunks overlap at " + length);
				assertEquals(QueryChunks.OVERLAP, c.end(i - 1) - c.start(i), "overlap at " + length);
			}
		}
	}
}
