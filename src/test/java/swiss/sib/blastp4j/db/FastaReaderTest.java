package swiss.sib.blastp4j.db;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

import org.junit.jupiter.api.Test;

class FastaReaderTest {

	static SequenceDatabase parse(String fasta) throws IOException {
		byte[] bytes = fasta.getBytes(StandardCharsets.US_ASCII);
		return FastaReader.read(new ByteArrayInputStream(bytes), bytes.length);
	}

	@Test
	void layoutMatchesPsq() throws IOException {
		SequenceDatabase db = parse(">sp|P1|A one\nMAC\nDE\n>sp|P2|B\nKW\n");
		assertEquals(2, db.size());
		assertEquals(7, db.totalResidues());
		assertEquals("sp|P1|A", db.id(0));
		assertEquals("sp|P2|B", db.id(1));
		assertEquals("MACDE", db.letters(0));
		assertEquals("KW", db.letters(1));
		// sentinel, M A C D E, sentinel, K W, sentinel (NCBIstdaa codes)
		assertArrayEquals(new byte[] { 0, 12, 1, 3, 4, 5, 0, 10, 20, 0 }, db.residues());
		assertEquals(1, db.start(0));
		assertEquals(7, db.start(1));
	}

	@Test
	void lowercaseCrlfUnknownLettersEmptySequenceAndNoFinalNewline() throws IOException {
		SequenceDatabase db = parse(">a\r\nma c\r\nu1o\r\n>empty\r\n>b x\r\nbzj*");
		assertEquals(3, db.size());
		assertEquals("a", db.id(0));
		assertEquals("MACUXO", db.letters(0)); // '1' is not a letter: X
		assertEquals(0, db.length(1));
		assertEquals("BZJ*", db.letters(2));
		assertEquals(10, db.totalResidues());
	}

	@Test
	void alphabetIsNcbistdaa() {
		assertEquals(28, Alphabet.SIZE);
		assertEquals(1, Alphabet.encode('A'));
		assertEquals(21, Alphabet.encode('X'));
		assertEquals(24, Alphabet.encode('U'));
		assertEquals(26, Alphabet.encode('O'));
		assertEquals(27, Alphabet.encode('J'));
		assertEquals(25, Alphabet.encode('*'));
		assertEquals('W', Alphabet.decode(20));
	}
}
