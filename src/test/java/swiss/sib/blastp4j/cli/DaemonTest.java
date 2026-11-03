package swiss.sib.blastp4j.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * The daemon: jobs posted over HTTP give blastp's output for their FASTA, their files move from input to output, and
 * the ids continue after a restart.
 */
class DaemonTest {
	@TempDir
	Path dir;

	private final HttpClient http = HttpClient.newHttpClient();

	Daemon start() throws IOException {
		return new Daemon(dir, Map.of("sample", Path.of(BlastpTest.DB)), new InetSocketAddress("127.0.0.1", 0), 2,
				1 << 20, new PrintStream(new ByteArrayOutputStream()));
	}

	HttpResponse<String> post(Daemon d, String params, String body) throws Exception {
		return http.send(HttpRequest.newBuilder(URI.create("http://127.0.0.1:" + d.port() + "/jobs?" + params))
				.POST(HttpRequest.BodyPublishers.ofString(body)).build(), HttpResponse.BodyHandlers.ofString());
	}

	HttpResponse<String> get(Daemon d, String path) throws Exception {
		return http.send(HttpRequest.newBuilder(URI.create("http://127.0.0.1:" + d.port() + path)).build(),
				HttpResponse.BodyHandlers.ofString());
	}

	/** The job's response once it is no longer queued or running. */
	HttpResponse<String> await(Daemon d, String id) throws Exception {
		for (int i = 0; i < 1200; i++) {
			HttpResponse<String> r = get(d, "/jobs/" + id);
			if (r.statusCode() != 202) {
				return r;
			}
			Thread.sleep(100);
		}
		throw new AssertionError("job " + id + " did not finish");
	}

	@Test
	void jobsAsBlastp() throws Exception {
		String queries = Files.readString(Path.of(BlastpTest.QUERIES));
		try (Daemon d = start()) {
			assertEquals("sample\t400\t", get(d, "/databases").body().substring(0, "sample\t400\t".length()));
			HttpResponse<String> first = post(d, "db=sample&outfmt=6+std+btop", queries);
			assertEquals(201, first.statusCode(), first.body());
			assertEquals("1", first.body().trim());
			assertEquals("/jobs/1", first.headers().firstValue("Location").orElseThrow());
			String second = post(d, "db=sample&outfmt=6%20std%20btop", queries).body().trim();
			assertEquals("2", second);
			HttpResponse<String> r = await(d, "1");
			assertEquals(200, r.statusCode(), r.body());
			assertEquals(BlastpTest.expected("sample_std_btop.tsv.gz"), r.body());
			assertEquals(BlastpTest.expected("sample_std_btop.tsv.gz"), await(d, second).body());
			assertTrue(Files.exists(dir.resolve("output/1.tsv")));
			assertFalse(Files.exists(dir.resolve("input/1.fasta")));
			assertFalse(Files.exists(dir.resolve("input/1.options")));
		}
		// The ids go on from the largest in the directories.
		try (Daemon d = start()) {
			assertEquals("3", post(d, "db=sample", ">q\nMKTAYIAKQRQISFVKSHFSRQ\n").body().trim());
			assertEquals(200, await(d, "3").statusCode());
		}
	}

	@Test
	void badJobs() throws Exception {
		try (Daemon d = start()) {
			assertEquals(400, post(d, "db=other", ">q\nMKTAYIAKQR\n").statusCode());
			assertEquals(400, post(d, "db=sample&gapopen=9", ">q\nMKTAYIAKQR\n").statusCode());
			assertEquals(400, post(d, "db=sample&num_threads=4", ">q\nMKTAYIAKQR\n").statusCode());
			assertEquals(400, post(d, "db=sample&outfmt=5", ">q\nMKTAYIAKQR\n").statusCode());
			assertEquals(400, post(d, "db=sample", "").statusCode());
			assertEquals(404, get(d, "/jobs/99").statusCode());
		}
		try (var files = Files.list(dir.resolve("input"))) {
			assertEquals(0, files.count());
		}
	}

	@Test
	void queuedJobsRunAfterARestart() throws Exception {
		Files.createDirectories(dir.resolve("input"));
		Files.createDirectories(dir.resolve("output"));
		Files.writeString(dir.resolve("input/7.fasta"), Files.readString(Path.of(BlastpTest.QUERIES)));
		Files.writeString(dir.resolve("input/7.options"), "db=sample\noutfmt=6 std btop\n", StandardCharsets.UTF_8);
		Files.writeString(dir.resolve("input/8.fasta"), ">partial\nMKT\n");
		Files.writeString(dir.resolve("output/5.tsv"), "");
		try (Daemon d = start()) {
			assertEquals(BlastpTest.expected("sample_std_btop.tsv.gz"), await(d, "7").body());
			assertEquals("9", post(d, "db=sample", ">q\nMKTAYIAKQRQISFVKSHFSRQ\n").body().trim());
			assertFalse(Files.exists(dir.resolve("input/8.fasta")));
		}
	}
}
