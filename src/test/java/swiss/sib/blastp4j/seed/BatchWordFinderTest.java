package swiss.sib.blastp4j.seed;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * The batch scan's hidden copies ({@link BatchWordFinder#of}) against {@link BatchWordFinder} itself: the copy for a
 * table with narrow (char) entries, and the copy for a batch too long for them, which reads the int entries.
 */
class BatchWordFinderTest {
	@Test
	void theCharCopyFindsWhatTheIntScanFinds() {
		Random random = new Random(3);
		byte[][] queries = queries(random, 24, 300);
		BatchWordFinder.Table table = table(queries);
		assertNotNull(table.narrow());
		compare(random, queries, table);
	}

	@Test
	void theCharCopyReadsTheNarrowEntries() {
		// A table whose narrow entries are all 0: if the copy reads them, it finds other alignments than the int scan.
		Random random = new Random(5);
		byte[][] queries = queries(random, 24, 300);
		BatchWordFinder.Table t = table(queries);
		BatchWordFinder.Table wrong = new BatchWordFinder.Table(t.start(), t.entry(), new char[t.narrow().length],
				t.query(), t.base(), t.size(), t.cells());
		LookupTable[] tables = tables(queries);
		SeedParameters[] parameters = parameters(queries.length);
		BatchScan copy = BatchWordFinder.of(queries, tables, parameters, wrong);
		BatchScan plain = new BatchWordFinder(queries, tables, parameters, wrong);
		boolean differ = false;
		for (int round = 0; round < 50 && !differ; round++) {
			byte[] q = queries[random.nextInt(queries.length)];
			byte[] subject = new byte[q.length + 2];
			System.arraycopy(q, 0, subject, 1, q.length);
			List<UngappedHsp>[] a = lists(queries.length), b = lists(queries.length);
			copy.find(subject, 1, q.length, a);
			plain.find(subject, 1, q.length, b);
			differ = !java.util.Arrays.equals(a, b);
		}
		assertTrue(differ, "the copy did not read the narrow entries");
	}

	@Test
	void aBatchTooLongForCharsIsScannedWithInts() {
		Random random = new Random(4);
		// Over 65,535 letters in all: the concatenation does not fit a char.
		byte[][] queries = queries(random, 3, 30_000);
		BatchWordFinder.Table table = table(queries);
		assertNull(table.narrow());
		compare(random, queries, table);
	}

	@Test
	void theIntDiagonalsOfALongSubjectFindWhatThe16BitOnesFind() {
		Random random = new Random(6);
		byte[][] queries = queries(random, 24, 300);
		BatchWordFinder.Table table = table(queries);
		LookupTable[] tables = tables(queries);
		SeedParameters[] parameters = parameters(queries.length);
		BatchWordFinder narrow = new BatchWordFinder(queries, tables, parameters, table);
		BatchWordFinder wide = new BatchWordFinder(queries, tables, parameters, table);
		int found = 0;
		// Enough subjects for the 16-bit offset to start again several times; the last one too long for 16 bits.
		for (int round = 0; round <= 300; round++) {
			int length = round < 300 ? 50 + random.nextInt(400) : 40_000;
			byte[] subject = new byte[length + 2];
			for (int at = 0; at < length;) {
				byte[] q = queries[random.nextInt(queries.length)];
				int from = random.nextInt(q.length), n = Math.min(length - at, 1 + random.nextInt(q.length - from));
				for (int i = 0; i < n; i++) {
					subject[1 + at + i] = random.nextInt(4) > 0 ? q[from + i] : (byte) (1 + random.nextInt(20));
				}
				at += n;
			}
			List<UngappedHsp>[] a = lists(queries.length), b = lists(queries.length);
			narrow.find(subject, 1, length, a);
			wide.findLong(subject, 1, length, b);
			for (int i = 0; i < queries.length; i++) {
				assertEquals(b[i], a[i], "round " + round + ", query " + i);
				found += a[i].size();
			}
		}
		assertTrue(found > 100, "alignments found: " + found);
	}

	/** The copy is a hidden class, and its alignments, query by query and in order, are those of the int scan. */
	private static void compare(Random random, byte[][] queries, BatchWordFinder.Table table) {
		LookupTable[] tables = tables(queries);
		SeedParameters[] parameters = parameters(queries.length);
		BatchScan copy = BatchWordFinder.of(queries, tables, parameters, table);
		BatchScan plain = new BatchWordFinder(queries, tables, parameters, table);
		assertTrue(copy.getClass().isHidden(), copy.getClass().getName());
		assertFalse(copy.getClass() == plain.getClass());
		int found = 0;
		for (int round = 0; round < 200; round++) {
			// A subject made of pieces of the queries with mutations, so that there are hits to extend, between
			// sentinels as in the database.
			byte[] q = queries[random.nextInt(queries.length)];
			int length = 50 + random.nextInt(400);
			byte[] subject = new byte[length + 2];
			int from = random.nextInt(Math.max(1, q.length - length));
			for (int i = 0; i < length; i++) {
				subject[i + 1] = from + i < q.length && random.nextInt(4) > 0 ? q[from + i]
						: (byte) (1 + random.nextInt(20));
			}
			List<UngappedHsp>[] a = lists(queries.length), b = lists(queries.length);
			copy.find(subject, 1, length, a);
			plain.find(subject, 1, length, b);
			for (int i = 0; i < queries.length; i++) {
				assertEquals(b[i], a[i], "round " + round + ", query " + i);
				found += a[i].size();
			}
		}
		assertTrue(found > 100, "alignments found: " + found);
	}

	private static byte[][] queries(Random random, int n, int length) {
		byte[][] queries = new byte[n][];
		for (int i = 0; i < n; i++) {
			queries[i] = new byte[length / 2 + random.nextInt(length)];
			for (int j = 0; j < queries[i].length; j++) {
				queries[i][j] = (byte) (1 + random.nextInt(20));
			}
		}
		return queries;
	}

	private static LookupTable[] tables(byte[][] queries) {
		LookupTable[] tables = new LookupTable[queries.length];
		for (int i = 0; i < queries.length; i++) {
			tables[i] = LookupTable.of(queries[i], 11);
		}
		return tables;
	}

	private static BatchWordFinder.Table table(byte[][] queries) {
		return BatchWordFinder.Table.of(queries, tables(queries), 40);
	}

	/** BLAST's window (40) and threshold (11), with a low X-drop and cutoff so that there are many alignments. */
	private static SeedParameters[] parameters(int n) {
		SeedParameters[] p = new SeedParameters[n];
		for (int i = 0; i < n; i++) {
			p[i] = new SeedParameters(11, 40, 16, 38, 20);
		}
		return p;
	}

	@SuppressWarnings("unchecked")
	private static List<UngappedHsp>[] lists(int n) {
		List<UngappedHsp>[] lists = new List[n];
		for (int i = 0; i < n; i++) {
			lists[i] = new ArrayList<>();
		}
		return lists;
	}
}
