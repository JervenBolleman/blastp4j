package swiss.sib.blastp4j.db;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.ReferenceData;

/** Queries are read with the same reader; their ids and lengths are BLAST's qseqid and qlen. */
class QueryReadingTest {

	@Test
	void queryIdsAndLengthsMatchTheReference() throws IOException {
		SequenceDatabase queries;
		try (InputStream in = QueryReadingTest.class.getResourceAsStream("/queries/uniprot_human.fasta")) {
			queries = FastaReader.read(in, 1 << 16);
		}
		Map<String, Integer> ours = new TreeMap<>();
		for (int i = 0; i < queries.size(); i++) {
			ours.put(queries.id(i), queries.length(i));
		}
		Map<String, Integer> blast = new TreeMap<>();
		for (String[] r : ReferenceData.load("human_vs_sprot.cbs0.tsv.gz")) {
			blast.put(r[ReferenceData.QSEQID], Integer.parseInt(r[ReferenceData.QLEN]));
		}
		assertEquals(24, ours.size());
		assertEquals(blast, ours);
	}
}
