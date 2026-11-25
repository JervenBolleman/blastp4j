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
	@Test
	void aSplitQueryGivesBlastpsHits() throws IOException, InterruptedException {
		final SequenceDatabase db = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
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
}
