package swiss.sib.blastp4j.stats;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;

class UngappedKarlinTest {
	/** blastp prints the parameters with three significant digits ({@code %#8.3g}). */
	static String blast(double v) {
		return String.format(Locale.ROOT, "%.3g", v);
	}

	@Test
	void sameAsBlastForEveryReferenceQuery() throws IOException {
		Map<String, String[]> expected = new HashMap<>();
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				getClass().getResourceAsStream("/reference/ungapped_karlin.tsv"), StandardCharsets.US_ASCII))) {
			in.readLine();
			for (String line; (line = in.readLine()) != null;) {
				String[] f = line.split("\t");
				expected.put(f[0], f);
			}
		}
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		assertEquals(queries.size(), expected.size());
		List<String> wrong = new ArrayList<>();
		for (int q = 0; q < queries.size(); q++) {
			UngappedKarlin k = UngappedKarlin.of(queries.sequence(q));
			String[] e = expected.get(queries.id(q));
			String ours = blast(k.lambda()) + " " + blast(k.k()) + " " + blast(k.h());
			String theirs = e[1] + " " + e[2] + " " + e[3];
			if (!ours.equals(theirs)) {
				wrong.add(queries.id(q) + ": " + ours + ", blastp " + theirs);
			}
		}
		assertTrue(wrong.isEmpty(), String.join("\n", wrong));
	}

	/** The classic BLOSUM62 values for the Robinson frequencies on both sides: lambda 0.3176, K 0.134, H 0.401. */
	@Test
	void backgroundAgainstItself() {
		StringBuilder sb = new StringBuilder();
		String letters = "ACDEFGHIKLMNPQRSTVWY";
		double[] counts = { 78.05, 19.25, 53.64, 62.95, 38.56, 73.77, 21.99, 51.42, 57.44, 90.19, 22.43, 44.87, 52.03,
				42.64, 51.29, 71.20, 58.41, 64.41, 13.30, 32.16 };
		for (int i = 0; i < letters.length(); i++) {
			sb.repeat(letters.charAt(i), (int) Math.round(counts[i] * 100));
		}
		UngappedKarlin k = UngappedKarlin.of(swiss.sib.blastp4j.db.Alphabet.encode(sb));
		assertEquals(0.3176, k.lambda(), 0.0001);
		assertEquals(0.134, k.k(), 0.001);
		assertEquals(0.401, k.h(), 0.001);
	}
}
