package swiss.sib.blastp4j.search;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.output.TabularFormat;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * A query long enough for BLAST to split it into chunks gives blastp's hits. Titin (34,350 residues) against the 24
 * human proteins is a small case of {@link QueryChunks}: of the five HSPs blastp reports, the one at query
 * 33,034-33,244 is there only because the third chunk's cutoff score is lower than the whole query's (31 against 34),
 * so a search of the query in one piece reports four.
 * <p>
 * Reference: {@code blastp -query titin.fasta -db uniprot_human.fasta -num_threads 1} with BLAST+ 2.17.0 and the
 * columns of {@code reference/README.md}.
 */
class QuerySplitReferenceTest {
	private static final String HUMAN = "src/test/resources/queries/uniprot_human.fasta";

	@Test
	void aSplitQueryGivesBlastpsHits() throws IOException, InterruptedException {
		final SequenceDatabase db = FastaReader.read(Path.of(HUMAN));
		final SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/titin.fasta"));
		final List<String> expected = Files.readAllLines(Path.of("src/test/resources/reference/titin_vs_human.tsv"));
		final byte[] query = queries.sequence(0);
		assertEquals(34350, query.length);
		assertEquals(3, QueryChunks.of(query.length).count(), "chunks of the query");

		final Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		final List<String> ours = new ArrayList<>();
		final int threads = Runtime.getRuntime().availableProcessors();
		try (ExecutorService pool = Executors.newFixedThreadPool(threads)) {
			final BlastpSearch search = new BlastpSearch(db, statistics, 10, true);
			for (BlastpSearch.Hit h : search.search(query, pool, threads)) {
				for (Hsp hsp : h.hsps()) {
					ours.add(TabularFormat.row(queries.id(0), query, db, h.subject(), hsp, statistics));
				}
			}
		}
		assertEquals(String.join("\n", expected), String.join("\n", ours));
	}

	/**
	 * The same query sharing a batch, which BLAST splits over the batch's queries together: SODC_HUMAN (154 residues)
	 * and titin are one batch of 34,504 residues, cut into three chunks of 11,568. SODC lies within the first, so
	 * that chunk holds the whole of it and a 11,414-residue piece of titin, and the other two hold pieces of titin
	 * that start 11,314 and 22,682 residues into it.
	 */
	@Test
	void aSplitBatchGivesBlastpsHits() throws IOException, InterruptedException {
		final SequenceDatabase db = FastaReader.read(Path.of(HUMAN));
		final SequenceDatabase human = FastaReader.read(Path.of(HUMAN));
		final SequenceDatabase titin = FastaReader.read(Path.of("src/test/resources/queries/titin.fasta"));
		final List<String> expected = Files
				.readAllLines(Path.of("src/test/resources/reference/sodc_titin_vs_human.tsv"));
		int sodc = 0;
		while (!human.id(sodc).equals("sp|P00441|SODC_HUMAN")) {
			sodc++;
		}
		final List<byte[]> batch = List.of(human.sequence(sodc), titin.sequence(0));
		final List<String> ids = List.of(human.id(sodc), titin.id(0));
		assertEquals(3, QueryChunks.of(154 + 34350).count(), "chunks of the batch");

		final Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		final List<String> ours = new ArrayList<>();
		final int threads = Runtime.getRuntime().availableProcessors();
		try (ExecutorService pool = Executors.newFixedThreadPool(threads)) {
			final BlastpSearch search = new BlastpSearch(db, statistics, 10, true);
			final QueryBatches batches = QueryBatches
					.blastp(new int[] { batch.get(0).length, batch.get(1).length });
			final List<List<BlastpSearch.Hit>> prelim = search.prelim(batch, pool, threads);
			for (int q = 0; q < batch.size(); q++) {
				final boolean earlyTermination = batches.next();
				final BlastpSearch.Result result = search.finish(batch.get(q), prelim.get(q), pool, threads,
						earlyTermination);
				batches.done(result.filledToCutoff());
				for (BlastpSearch.Hit h : result.hits()) {
					for (Hsp hsp : h.hsps()) {
						ours.add(TabularFormat.row(ids.get(q), batch.get(q), db, h.subject(), hsp, statistics));
					}
				}
			}
		}
		assertEquals(String.join("\n", expected), String.join("\n", ours));
	}
}
