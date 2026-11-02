package swiss.sib.blastp4j.search;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/** Where BLAST's composition-based stage may skip subjects, given its query batches. */
class QueryBatchesTest {
	@Test
	void onlyTheLastQueryOfABatchAfterFilledOnes() {
		// Batches of 100 residues: [60, 50], [120], [30, 30, 30] (the input ends before 100).
		QueryBatches b = new QueryBatches(new int[] { 60, 50, 120, 30, 30, 30 }, 100);
		assertFalse(b.next());
		b.done(true);
		assertTrue(b.next(), "last of the first batch, the earlier query filled");
		b.done(false);
		assertTrue(b.next(), "alone in its batch");
		b.done(false);
		assertFalse(b.next());
		b.done(true);
		assertFalse(b.next());
		b.done(false);
		assertFalse(b.next(), "last of the input, but an earlier query of its batch did not fill");
		b.done(true);
	}

	@Test
	void aLongQueryIsAloneInItsBatch() {
		QueryBatches b = new QueryBatches(new int[] { 20000, 10 }, QueryBatches.BLASTP_BATCH_SIZE);
		assertTrue(b.next());
		b.done(false);
		assertTrue(b.next(), "alone in the last batch");
		b.done(false);
	}
}
