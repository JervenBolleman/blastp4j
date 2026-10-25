package swiss.sib.blastp4j.search;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import swiss.sib.blastp4j.ReferenceData;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.output.TabularFormat;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * The complete search reproduces blastp's tabular output line for line, in the same order, for three of the reference
 * queries, without composition-based statistics and with blastp's default mode 2 (all 24 are compared by
 * {@link GappedBaseline}). Skipped when Swiss-Prot is not in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 * <p>
 * CRYAA has subjects that need the restricted alignment and its exact redo; UBC (four ubiquitin repeats) has many
 * HSPs per subject, where the containment rules matter.
 */
class BlastpSearchReferenceTest {
	static SequenceDatabase db, queries;

	@BeforeAll
	static void read() throws IOException {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		assumeTrue(Files.exists(dir.resolve("uniprot_sprot.fasta")), "Swiss-Prot not found in " + dir);
		db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
	}

	@ParameterizedTest
	@CsvSource({ "sp|P00441|SODC_HUMAN,cbs0", "sp|P02489|CRYAA_HUMAN,cbs0", "sp|P0CG48|UBC_HUMAN,cbs0",
			"sp|P00441|SODC_HUMAN,cbs2", "sp|P02489|CRYAA_HUMAN,cbs2", "sp|P0CG48|UBC_HUMAN,cbs2" })
	void sameLinesAsBlastp(String qid, String mode) throws Exception {
		List<String> expected = new ArrayList<>();
		for (String[] r : ReferenceData.load("human_vs_sprot." + mode + ".tsv.gz")) {
			if (r[0].equals(qid)) {
				expected.add(String.join("\t", r));
			}
		}
		int q = 0;
		while (!queries.id(q).equals(qid)) {
			q++;
		}
		byte[] query = queries.sequence(q);
		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		List<String> ours = new ArrayList<>();
		int threads = Runtime.getRuntime().availableProcessors();
		try (ExecutorService pool = Executors.newFixedThreadPool(threads)) {
			BlastpSearch search = new BlastpSearch(db, statistics, 10, mode.equals("cbs2"));
			for (BlastpSearch.Hit h : search.search(query, pool, threads)) {
				for (Hsp hsp : h.hsps()) {
					ours.add(TabularFormat.row(qid, db.id(h.subject()), query, db.residues(), db.start(h.subject()),
							db.length(h.subject()), hsp, statistics));
				}
			}
		}
		assertEquals(String.join("\n", expected), String.join("\n", ours));
	}
}
