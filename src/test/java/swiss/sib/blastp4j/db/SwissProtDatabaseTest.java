package swiss.sib.blastp4j.db;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.ReferenceData;

/**
 * Reads UniProtKB/Swiss-Prot 2026_03 and compares it with the BLAST database built from the same file (see
 * reference/README.md). Skipped when the data is not present; the directory is {@code ~/data/uniprot} or the system
 * property {@code blastp4j.sprot}.
 */
class SwissProtDatabaseTest {
	static final Path DIR = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
	static final Path FASTA = DIR.resolve("uniprot_sprot.fasta");
	static final Path PSQ = DIR.resolve("blastdb/sprot.psq");

	static SequenceDatabase db;

	@BeforeAll
	static void read() throws IOException {
		assumeTrue(Files.exists(FASTA) && Files.exists(PSQ), "Swiss-Prot not found in " + DIR);
		db = FastaReader.read(FASTA);
	}

	@Test
	void sizesMatchMakeblastdb() {
		assertEquals(ReferenceData.SPROT_SEQUENCES, db.size());
		assertEquals(ReferenceData.SPROT_RESIDUES, db.totalResidues());
	}

	@Test
	void residuesAreIdenticalToThePsqFile() throws IOException {
		byte[] psq = Files.readAllBytes(PSQ);
		byte[] ours = db.residues();
		assertEquals(psq.length, ours.length, "length");
		int first = java.util.Arrays.mismatch(psq, ours);
		assertEquals(-1, first, "first differing byte");
	}

	@Test
	void everyReferenceSubjectHasItsIdAndLength() throws IOException {
		Map<String, Integer> lengths = new HashMap<>();
		for (int i = 0; i < db.size(); i++) {
			lengths.put(db.id(i), db.length(i));
		}
		assertEquals(db.size(), lengths.size(), "ids are unique");
		List<String> wrong = new ArrayList<>();
		for (String name : new String[] { "human_vs_sprot.cbs0.tsv.gz", "human_vs_sprot.cbs2.tsv.gz" }) {
			for (String[] r : ReferenceData.load(name)) {
				Integer length = lengths.get(r[ReferenceData.SSEQID]);
				if (length == null || length != Integer.parseInt(r[ReferenceData.SLEN])) {
					wrong.add(r[ReferenceData.SSEQID] + " length " + length + ", expected " + r[ReferenceData.SLEN]);
				}
			}
		}
		assertTrue(wrong.isEmpty(), wrong.size() + " subjects differ, e.g. " + wrong.subList(0, Math.min(5, wrong.size())));
	}
}
