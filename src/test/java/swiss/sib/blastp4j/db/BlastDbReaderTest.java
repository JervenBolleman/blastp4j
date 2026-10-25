package swiss.sib.blastp4j.db;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * Reading BLAST databases made by makeblastdb 2.17.0 (see src/test/resources/blastdb/README.md): the sequences must be
 * those of the FASTA file they were made from, and the identifiers those blastp prints.
 */
class BlastDbReaderTest {
	static final Path DIR = Path.of("src/test/resources/blastdb");
	static final Path SAMPLE = Path.of("src/test/resources/seg/sprot_sample.fasta.gz");

	static void assertSameSequences(SequenceDatabase expected, SequenceDatabase actual) {
		assertEquals(expected.size(), actual.size());
		assertEquals(expected.totalResidues(), actual.totalResidues());
		for (int i = 0; i < expected.size(); i++) {
			assertArrayEquals(expected.sequence(i), actual.sequence(i), "sequence " + i);
		}
		assertEquals(expected.maxLength(), actual.maxLength());
	}

	@Test
	void identifiersWithParseSeqids() throws IOException {
		SequenceDatabase fasta = FastaReader.read(DIR.resolve("ids.fa"));
		SequenceDatabase db = BlastDbReader.read(DIR.resolve("ids_p"));
		assertSameSequences(fasta, db);
		// Identifiers as blastp prints them; the FASTA reader parses headers as makeblastdb -parse_seqids does.
		for (int i = 0; i < db.size(); i++) {
			assertEquals(fasta.id(i), db.id(i));
			assertEquals(fasta.accession(i, false), db.accession(i, false));
			assertEquals(fasta.accession(i, true), db.accession(i, true));
			assertEquals(fasta.title(i), db.title(i));
		}
		assertEquals("tr|A0A011|A0A011_HUMAN", db.id(1));
		assertEquals("pdb|1ABC|A", db.id(4));
		assertEquals("1ABC_A", db.accession(4, true));
	}

	@Test
	void identifiersWithoutParseSeqids() throws IOException {
		SequenceDatabase db = BlastDbReader.read(DIR.resolve("ids_n"));
		// The identifier is gnl|BL_ORD_ID|oid; blastp prints the first word of the definition line instead.
		assertEquals("lcl|myseq", db.id(5));
		assertEquals("lcl|myseq", db.accession(5, true));
		assertEquals("lcl|myseq local id", db.title(5));
	}

	@Test
	void multipleVolumes() throws IOException {
		SequenceDatabase fasta = FastaReader.read(SAMPLE);
		SequenceDatabase db = BlastDbReader.read(DIR.resolve("multi"));
		assertSameSequences(fasta, db);
		for (int i = 0; i < db.size(); i++) {
			assertEquals(fasta.id(i), db.id(i));
			assertEquals(fasta.accession(i, true), db.accession(i, true));
			assertEquals(fasta.title(i), db.title(i));
		}
		assertEquals(6, db.volumes());
	}

	@Test
	void exists(@TempDir Path tmp) throws IOException {
		assertTrue(BlastDbReader.exists(DIR.resolve("ids_p")));
		assertTrue(BlastDbReader.exists(DIR.resolve("multi")));
		assertFalse(BlastDbReader.exists(DIR.resolve("ids.fa")));
		Files.writeString(tmp.resolve("restricted.pal"), "DBLIST " + DIR.toAbsolutePath().resolve("ids_p")
				+ "\nOIDLIST some.msk\n");
		assertThrows(IOException.class, () -> BlastDbReader.read(tmp.resolve("restricted")));
		assertThrows(IOException.class, () -> BlastDbReader.read(tmp.resolve("missing")));
	}

	@Test
	void swissProt() throws IOException {
		Path db = SwissProtDatabaseTest.DIR.resolve("blastdb/sprot");
		Path fasta = SwissProtDatabaseTest.FASTA;
		Assumptions.assumeTrue(BlastDbReader.exists(db) && Files.exists(fasta), "Swiss-Prot not downloaded");
		SequenceDatabase expected = FastaReader.read(fasta);
		SequenceDatabase actual = BlastDbReader.read(db);
		assertSameSequences(expected, actual);
		for (int i = 0; i < actual.size(); i += 97) {
			assertEquals(expected.id(i), actual.id(i));
			assertEquals(expected.accession(i, true), actual.accession(i, true));
			assertEquals(expected.title(i), actual.title(i));
		}
	}
}
