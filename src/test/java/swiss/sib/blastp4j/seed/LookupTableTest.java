package swiss.sib.blastp4j.seed;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.blastp4j.stats.Blosum62;

class LookupTableTest {
	/** Every word of the alphabet against every query word, without pruning. */
	static int[][] bruteForce(byte[] query, int threshold) {
		int[][] cells = new int[LookupTable.WORDS][];
		for (int x = 0; x < Blosum62.SIZE; x++) {
			for (int y = 0; y < Blosum62.SIZE; y++) {
				for (int z = 0; z < Blosum62.SIZE; z++) {
					List<Integer> hits = new ArrayList<>();
					for (int q = 0; q + 3 <= query.length; q++) {
						int score = Blosum62.score(query[q], x) + Blosum62.score(query[q + 1], y)
								+ Blosum62.score(query[q + 2], z);
						boolean exact = query[q] == x && query[q + 1] == y && query[q + 2] == z;
						if (score >= threshold || exact) {
							hits.add(q);
						}
					}
					cells[(x << 10) | (y << 5) | z] = hits.stream().mapToInt(Integer::intValue).toArray();
				}
			}
		}
		return cells;
	}

	@Test
	void sameAsBruteForce() {
		// Includes a repeated word, low-scoring words (AAA scores 12, GSA only 10: exact match only) and X.
		byte[] query = Alphabet.encode("MKWVTFISLLFLFSSAYSGSAAAXGSAWWWAAA");
		LookupTable table = LookupTable.of(query, 11);
		int[][] expected = bruteForce(query, 11);
		int total = 0;
		for (int w = 0; w < LookupTable.WORDS; w++) {
			int[] ours = table.offsets(w);
			Arrays.sort(ours);
			int[] theirs = expected[w] == null ? new int[0] : expected[w];
			assertArrayEquals(theirs, ours, "word " + w);
			total += theirs.length;
		}
		assertEquals(total, table.size());
	}
}
