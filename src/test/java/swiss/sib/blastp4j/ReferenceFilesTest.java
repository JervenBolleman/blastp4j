package swiss.sib.blastp4j;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

/** The BLAST+ reference results are present and have the expected shape (see reference/README.md). */
class ReferenceFilesTest {

	@Test
	void referencesHaveAllColumns() throws IOException {
		for (String name : new String[] { "human_vs_sprot.cbs2.tsv.gz", "human_vs_sprot.cbs0.tsv.gz" }) {
			List<String[]> rows = ReferenceData.load(name);
			for (String[] row : rows) {
				assertEquals(16, row.length, name + ": " + String.join("\t", row));
			}
			assertEquals(24, rows.stream().map(r -> r[0]).distinct().count(), name + ": queries");
		}
		assertEquals(11_263, ReferenceData.load("human_vs_sprot.cbs2.tsv.gz").size());
		assertEquals(15_107, ReferenceData.load("human_vs_sprot.cbs0.tsv.gz").size());
	}
}
