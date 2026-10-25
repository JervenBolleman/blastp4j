package swiss.sib.blastp4j.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import swiss.sib.blastp4j.db.SeqId;
import swiss.sib.blastp4j.output.TabularFormat;

/**
 * The command line: option handling, and output identical to blastp's for the 24 reference queries against a small
 * database, the 400 Swiss-Prot sequences of the SEG sample (blastp run against {@code makeblastdb -parse_seqids} of
 * the same sequences; see reference/README.md).
 */
class BlastpTest {
	static final String QUERIES = "src/test/resources/queries/uniprot_human.fasta";
	static final String DB = "src/test/resources/seg/sprot_sample.fasta.gz";

	@TempDir
	Path tmp;

	static String expected(String name) throws IOException {
		try (InputStream in = new GZIPInputStream(BlastpTest.class.getResourceAsStream("/cli/" + name))) {
			return new String(in.readAllBytes(), StandardCharsets.US_ASCII);
		}
	}

	String run(String... args) throws IOException {
		Path out = tmp.resolve("out.tsv");
		String[] all = new String[args.length + 2];
		System.arraycopy(args, 0, all, 0, args.length);
		all[args.length] = "-out";
		all[args.length + 1] = out.toString();
		ByteArrayOutputStream err = new ByteArrayOutputStream();
		int status = Blastp.run(all, new PrintStream(new ByteArrayOutputStream()), new PrintStream(err));
		assertEquals(0, status, err.toString());
		return Files.readString(out, StandardCharsets.US_ASCII);
	}

	@Test
	void defaultsAsBlastp() throws IOException {
		assertEquals(expected("sample_std_btop.tsv.gz"),
				run("-query", QUERIES, "-db", DB, "-outfmt", "6 std btop", "-num_threads", "2"));
	}

	@Test
	void optionsAsBlastp() throws IOException {
		// -evalue above 10 also switches off the restricted alignment of the preliminary stage.
		assertEquals(expected("sample_cbs0_e1000.tsv.gz"),
				run("-query", QUERIES, "-db", DB, "-comp_based_stats", "0", "-evalue", "1000", "-max_target_seqs", "20",
						"-outfmt", "6 qseqid sseqid qcovhsp ppos positive gaps nident evalue bitscore score"));
	}

	@Test
	void blastDatabase() throws IOException {
		// The sample as a BLAST database of six volumes gives the same result as the FASTA file.
		assertEquals(expected("sample_std_btop.tsv.gz"), run("-query", QUERIES, "-db", "src/test/resources/blastdb/multi",
				"-outfmt", "6 std btop", "-num_threads", "2"));
		// Identifiers of every kind, with and without makeblastdb -parse_seqids.
		for (String db : new String[] { "ids_p", "ids_n" }) {
			Path dir = Path.of("src/test/resources/blastdb");
			assertEquals(Files.readString(dir.resolve(db + ".tsv")),
					run("-query", dir.resolve("q.fa").toString(), "-db", dir.resolve(db).toString(), "-evalue", "1000",
							"-outfmt", "6 qseqid sseqid sacc saccver stitle slen evalue bitscore"), db);
		}
	}

	@Test
	void parsesOptions() throws Blastp.UsageException {
		Blastp.Options o = Blastp.parse(new String[] { "-query", "q.fa", "-db", "d.fa", "-evalue", "1e-5",
				"-max_target_seqs", "50", "-max_hsps", "2", "-comp_based_stats", "F", "-num_threads", "3", "-task",
				"blastp", "-seg", "no" });
		assertEquals(1e-5, o.evalue());
		assertEquals(50, o.maxTargetSeqs());
		assertEquals(2, o.maxHsps());
		assertEquals(false, o.compositionBased());
		assertEquals(3, o.threads());
		assertEquals(TabularFormat.STANDARD, o.columns());
	}

	@Test
	void rejectsWhatIsNotImplemented() {
		assertThrows(Blastp.UsageException.class,
				() -> Blastp.parse(new String[] { "-query", "q", "-db", "d", "-matrix", "PAM30" }));
		assertThrows(Blastp.UsageException.class,
				() -> Blastp.parse(new String[] { "-query", "q", "-db", "d", "-outfmt", "5" }));
		assertThrows(Blastp.UsageException.class,
				() -> Blastp.parse(new String[] { "-query", "q", "-db", "d", "-outfmt", "6 qseqid staxids" }));
		assertThrows(Blastp.UsageException.class,
				() -> Blastp.parse(new String[] { "-query", "q", "-db", "d", "-comp_based_stats", "1" }));
		assertThrows(Blastp.UsageException.class, () -> Blastp.parse(new String[] { "-query", "q" }));
		assertThrows(Blastp.UsageException.class,
				() -> Blastp.parse(new String[] { "-query", "q", "-db", "d", "-evalue", "ten" }));
	}

	@Test
	void helpAndErrors() {
		ByteArrayOutputStream out = new ByteArrayOutputStream(), err = new ByteArrayOutputStream();
		assertEquals(0, Blastp.run(new String[] { "-help" }, new PrintStream(out), new PrintStream(err)));
		assertTrue(out.toString().contains("-max_target_seqs"));
		assertEquals(1, Blastp.run(new String[] { "-query", "q" }, new PrintStream(out), new PrintStream(err)));
		assertTrue(err.toString().contains("-query and -db are required"), err.toString());
		assertEquals(2, Blastp.run(new String[] { "-query", "missing.fa", "-db", DB }, new PrintStream(out),
				new PrintStream(err)));
	}

	@Test
	void accessions() {
		assertEquals("P69905", acc("sp|P69905|HBA_HUMAN", true));
		assertEquals("NP_000509.1", acc("ref|NP_000509.1|", true));
		assertEquals("NP_000509", acc("ref|NP_000509.1|", false));
		assertEquals("NP_222222", acc("gi|123456|ref|NP_222222.2|", false));
		assertEquals("1ABC_A", acc("pdb|1ABC|A", true));
		assertEquals("mydb:xyz", acc("gnl|mydb|xyz", true));
		assertEquals("myseq", acc("lcl|myseq", true));
		assertEquals("seq1", acc("seq1", true));
		assertEquals("myseq", SeqId.seqid(SeqId.parseFasta("lcl|myseq")));
		assertEquals("gi|123456|ref|NP_222222.2|", SeqId.seqid(SeqId.parseFasta("gi|123456|ref|NP_222222.2|")));
		assertEquals("tr|A0A011|A0A011_HUMAN", SeqId.seqid(SeqId.parseFasta("tr|A0A011|A0A011_HUMAN")));
	}

	static String acc(String id, boolean version) {
		return SeqId.accession(SeqId.parseFasta(id), version);
	}

}
