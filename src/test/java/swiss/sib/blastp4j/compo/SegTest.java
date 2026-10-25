package swiss.sib.blastp4j.compo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;

/**
 * SEG with the parameters of composition-based statistics masks exactly what BLAST+'s {@code segmasker -window 10
 * -locut 1.8 -hicut 2.1} masks, for 400 Swiss-Prot sequences (every 50th of the first 20,000; see
 * reference/README.md).
 */
class SegTest {
	static InputStream resource(String name) throws IOException {
		return new GZIPInputStream(SegTest.class.getResourceAsStream("/seg/" + name));
	}

	@Test
	void sameMasksAsSegmasker() throws IOException {
		SequenceDatabase db;
		try (InputStream in = resource("sprot_sample.fasta.gz")) {
			db = FastaReader.read(in, 200_000);
		}
		List<BitSet> expected = new ArrayList<>();
		try (BufferedReader in = new BufferedReader(
				new InputStreamReader(resource("sprot_sample.segmasker.txt.gz"), StandardCharsets.US_ASCII))) {
			BitSet current = null;
			for (String line; (line = in.readLine()) != null;) {
				if (line.startsWith(">")) {
					current = new BitSet();
					expected.add(current);
				} else {
					String[] f = line.split(" - ");
					current.set(Integer.parseInt(f[0].trim()), Integer.parseInt(f[1].trim()) + 1);
				}
			}
		}
		assertEquals(db.size(), expected.size());
		Seg seg = Seg.forCompositionBasedStatistics();
		int masked = 0;
		for (int i = 0; i < db.size(); i++) {
			BitSet ours = new BitSet();
			for (int[] s : seg.segments(db.withSentinels(i, null), 1, db.length(i))) {
				ours.set(s[0], s[1] + 1);
			}
			assertEquals(expected.get(i), ours, db.id(i));
			masked += ours.isEmpty() ? 0 : 1;
		}
		assertTrue(masked > 100, "the sample has masked sequences: " + masked);
	}
}
