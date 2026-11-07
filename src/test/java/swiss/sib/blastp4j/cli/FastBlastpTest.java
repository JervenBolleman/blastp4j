package swiss.sib.blastp4j.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * The fastblastp command: its output is BLAST+'s {@code blastp -task blastp-fast} for the 24 queries against the 400
 * sequences of the SEG sample (see reference/README.md).
 */
class FastBlastpTest {
	@TempDir
	Path tmp;

	String run(String... args) throws IOException {
		Path out = tmp.resolve("out.tsv");
		String[] all = new String[args.length + 2];
		System.arraycopy(args, 0, all, 0, args.length);
		all[args.length] = "-out";
		all[args.length + 1] = out.toString();
		ByteArrayOutputStream err = new ByteArrayOutputStream();
		int status = FastBlastp.run(all, new PrintStream(new ByteArrayOutputStream()), new PrintStream(err));
		assertEquals(0, status, err.toString());
		return Files.readString(out, StandardCharsets.US_ASCII);
	}

	@Test
	void sameAsBlastpFast() throws IOException {
		assertEquals(BlastpTest.expected("sample_fast_std_btop.tsv.gz"), run("-query", BlastpTest.QUERIES, "-db",
				BlastpTest.DB, "-outfmt", "6 std btop", "-num_threads", "2"));
	}

	@Test
	void sameAsBlastpFastWithoutCompositionBasedStatistics() throws IOException {
		assertEquals(BlastpTest.expected("sample_fast_cbs0_e1000.tsv.gz"),
				run("-query", BlastpTest.QUERIES, "-db", BlastpTest.DB, "-comp_based_stats", "0", "-evalue", "1000",
						"-max_target_seqs", "20", "-outfmt",
						"6 qseqid sseqid qcovhsp ppos positive gaps nident evalue bitscore score"));
	}

	@Test
	void acceptsOnlyBlastpFastsParameters() {
		ByteArrayOutputStream err = new ByteArrayOutputStream();
		assertEquals(1, FastBlastp.run(new String[] { "-query", "q", "-db", "d", "-word_size", "3" },
				new PrintStream(new ByteArrayOutputStream()), new PrintStream(err)), err.toString());
		assertTrue(err.toString().contains("-word_size 3 is not supported"), err.toString());
	}

	@Test
	void namesItselfInHelpAndVersion() {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		assertEquals(0, FastBlastp.run(new String[] { "-help" }, new PrintStream(out), System.err));
		assertTrue(out.toString().startsWith("Usage: fastblastp "), out.toString());
		out.reset();
		assertEquals(0, FastBlastp.run(new String[] { "-version" }, new PrintStream(out), System.err));
		assertTrue(out.toString().startsWith("fastblastp "), out.toString());
	}
}
