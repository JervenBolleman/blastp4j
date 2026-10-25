package swiss.sib.blastp4j;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;

/** The BLAST+ reference results are present and have the expected shape (see reference/README.md). */
class ReferenceFilesTest {

	static List<String[]> load(String name) throws IOException {
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(ReferenceFilesTest.class.getResourceAsStream("/reference/" + name)),
				StandardCharsets.US_ASCII))) {
			return in.lines().map(l -> l.split("\t", -1)).toList();
		}
	}

	@Test
	void referencesHaveAllColumns() throws IOException {
		for (String name : new String[] { "human_vs_sprot.cbs2.tsv.gz", "human_vs_sprot.cbs0.tsv.gz" }) {
			List<String[]> rows = load(name);
			for (String[] row : rows) {
				assertEquals(16, row.length, name + ": " + String.join("\t", row));
			}
			assertEquals(24, rows.stream().map(r -> r[0]).distinct().count(), name + ": queries");
		}
		assertEquals(11_263, load("human_vs_sprot.cbs2.tsv.gz").size());
		assertEquals(15_107, load("human_vs_sprot.cbs0.tsv.gz").size());
	}
}
