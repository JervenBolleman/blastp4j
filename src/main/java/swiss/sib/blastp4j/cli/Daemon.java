package swiss.sib.blastp4j.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.io.UncheckedIOException;
import java.io.Writer;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;

import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * {@code blastp4j daemon}: a search server. It opens its databases once, takes jobs over HTTP and runs them one at a
 * time, each on all threads.
 *
 * <pre>
 * blastp4j daemon -dir jobs -db sprot=/data/blastdb/sprot -db bean=proteome.fasta -port 8080
 * curl --data-binary @queries.fasta 'http://localhost:8080/jobs?db=sprot&amp;outfmt=6+std+btop'   # prints the job id
 * curl -i http://localhost:8080/jobs/1                                                    # 202 until done
 * </pre>
 *
 * A job is stored as {@code input/<id>.fasta} and {@code input/<id>.options} in the job directory. When it is done its
 * tabular output is {@code output/<id>.tsv} (or the error {@code output/<id>.error}) and its input is deleted. Job ids
 * count up from one more than the largest id in either directory, and the jobs left in {@code input} are run again at
 * startup. A job's queries are searched in blastp's batches, as {@code blastp4j -query} does with the job's FASTA file,
 * so its output is the one blastp gives for that file. Jobs are not merged into one batch: BLAST's composition-based
 * stage depends on the batch, so the output would no longer be blastp's for each job's file.
 */
public final class Daemon implements AutoCloseable {
	static final String USAGE = """
			Usage: blastp4j daemon -dir DIR -db NAME=DB [-db NAME=DB]... [options]

			Runs a search server: jobs are posted over HTTP and run one at a time on all threads.

			  -dir DIR           job directory: input/ holds the queued jobs, output/ the results
			  -db NAME=DB        a database jobs may name: a BLAST protein database or a FASTA file
			  -port N            HTTP port (default 8080; 0 for any free port)
			  -bind ADDRESS      address to listen on (default 127.0.0.1)
			  -num_threads N     threads for each job (default: number of processors)
			  -max_query_mb N    largest query FASTA accepted, in MB (default 64)

			HTTP:
			  GET  /databases          the databases: name, sequences, residues (tab-separated)
			  POST /jobs?db=NAME&...   a job: the body is the query FASTA; parameters as blastp's options
			                           (outfmt, evalue, max_target_seqs, max_hsps, comp_based_stats, and the fixed
			                           ones); answers 201 with the job id
			  GET  /jobs/ID            200 with the tabular output when done, 202 while queued or running,
			                           500 with the error when it failed, 404 for an unknown job
			""";

	/** The blastp options a job may set; the others ({@code -query}, {@code -out}, threads) are the daemon's. */
	private static final Set<String> JOB_OPTIONS;
	static {
		Set<String> options = new TreeSet<>(
				List.of("outfmt", "evalue", "max_target_seqs", "max_hsps", "comp_based_stats"));
		options.addAll(Blastp.FIXED.keySet());
		JOB_OPTIONS = Set.copyOf(options);
	}

	private static final Pattern JOB_FILE = Pattern.compile("(\\d+)\\..*");
	private static final Pattern JOB_PATH = Pattern.compile("/jobs/(\\d+)/?");

	/** An open database and its statistics. */
	private record Database(Path path, SequenceDatabase db, Statistics statistics) {
	}

	private enum State {
		QUEUED, RUNNING
	}

	private final Path input, output;
	private final Map<String, Database> databases;
	private final int threads;
	private final long maxQueryBytes;
	private final PrintStream log;
	private final AtomicLong nextId;
	private final Map<Long, State> pending = new ConcurrentHashMap<>();
	private final LinkedBlockingQueue<Long> queue = new LinkedBlockingQueue<>();
	private final ExecutorService pool;
	private final Thread worker;
	private final HttpServer server;

	/**
	 * Opens the databases, queues the jobs left in {@code dir/input} and starts the server.
	 *
	 * @param databases the database paths by the names jobs use
	 */
	public Daemon(Path dir, Map<String, Path> databases, InetSocketAddress address, int threads, long maxQueryBytes,
			PrintStream log) throws IOException {
		this.input = Files.createDirectories(dir.resolve("input"));
		this.output = Files.createDirectories(dir.resolve("output"));
		this.threads = threads;
		this.maxQueryBytes = maxQueryBytes;
		this.log = log;
		Map<String, Database> open = new LinkedHashMap<>();
		try {
			for (Map.Entry<String, Path> e : databases.entrySet()) {
				SequenceDatabase db = Blastp.openDatabase(e.getValue());
				open.put(e.getKey(), new Database(e.getValue(), db,
						new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues())));
				log.println("blastp4j daemon: database " + e.getKey() + ": " + e.getValue() + ", " + db.size()
						+ " sequences");
			}
		} catch (IOException | RuntimeException e) {
			for (Database d : open.values()) {
				d.db().close();
			}
			throw e;
		}
		this.databases = Map.copyOf(open);
		this.nextId = new AtomicLong(recover() + 1);
		this.pool = Executors.newFixedThreadPool(threads);
		this.worker = Thread.ofPlatform().name("blastp4j-jobs").start(this::work);
		this.server = HttpServer.create(address, 0);
		server.createContext("/databases", this::databases);
		server.createContext("/jobs", this::jobs);
		server.setExecutor(Executors.newVirtualThreadPerTaskExecutor());
		server.start();
		log.println("blastp4j daemon: listening on " + server.getAddress() + ", jobs in " + dir + ", next job "
				+ nextId.get());
	}

	/** The port the server listens on. */
	public int port() {
		return server.getAddress().getPort();
	}

	/**
	 * Queues the complete jobs left in the input directory, in id order, and deletes partial files. Returns the
	 * largest job id in the input and output directories, 0 if none.
	 */
	private long recover() throws IOException {
		long max = 0;
		TreeSet<Long> queued = new TreeSet<>(), inputs = new TreeSet<>();
		for (Path d : List.of(input, output)) {
			try (DirectoryStream<Path> files = Files.newDirectoryStream(d)) {
				for (Path f : files) {
					String name = f.getFileName().toString();
					Matcher m = JOB_FILE.matcher(name);
					if (!m.matches()) {
						continue;
					}
					long id = Long.parseLong(m.group(1));
					max = Math.max(max, id);
					if (name.endsWith(".part")) {
						Files.delete(f);
					} else if (d == input) {
						(name.endsWith(".options") ? queued : inputs).add(id);
					}
				}
			}
		}
		for (long id : queued) {
			if (Files.exists(fasta(id))) {
				pending.put(id, State.QUEUED);
				queue.add(id);
			} else {
				Files.deleteIfExists(options(id));
			}
		}
		// A query without its options file is from a submission that did not complete.
		inputs.removeAll(queued);
		for (long id : inputs) {
			Files.deleteIfExists(fasta(id));
		}
		if (!queued.isEmpty()) {
			log.println("blastp4j daemon: " + queue.size() + " jobs left from before queued again");
		}
		return max;
	}

	private Path fasta(long id) {
		return input.resolve(id + ".fasta");
	}

	private Path options(long id) {
		return input.resolve(id + ".options");
	}

	private Path result(long id) {
		return output.resolve(id + ".tsv");
	}

	private Path error(long id) {
		return output.resolve(id + ".error");
	}

	/** The worker: runs the queued jobs in turn. */
	private void work() {
		try {
			while (true) {
				long id = queue.take();
				pending.put(id, State.RUNNING);
				try {
					run(id);
				} finally {
					pending.remove(id);
				}
			}
		} catch (InterruptedException e) {
			// Closed; the jobs still in the input directory run at the next start.
		}
	}

	/** Runs job {@code id}: writes its output (or error) to the output directory, then deletes its input. */
	private void run(long id) throws InterruptedException {
		long start = System.nanoTime();
		Path part = output.resolve(id + ".tsv.part");
		try {
			Properties p = new Properties();
			try (Reader r = Files.newBufferedReader(options(id), StandardCharsets.UTF_8)) {
				p.load(r);
			}
			Database d = databases.get(p.getProperty("db"));
			if (d == null) {
				throw new IOException("database " + p.getProperty("db") + " is no longer served");
			}
			Map<String, String> params = new LinkedHashMap<>();
			for (String name : p.stringPropertyNames()) {
				params.put(name, p.getProperty(name));
			}
			Blastp.Options o = options(fasta(id), d, params);
			int size;
			try (SequenceDatabase queries = FastaReader.read(fasta(id));
					Writer out = Files.newBufferedWriter(part, StandardCharsets.US_ASCII)) {
				size = queries.size();
				Blastp.search(o, queries, d.db(), d.statistics(), pool, out, new PrefixStream(log, "job " + id + ": "));
			}
			Files.move(part, result(id), StandardCopyOption.ATOMIC_MOVE);
			log.printf(Locale.ROOT, "blastp4j daemon: job %d done, %d queries against %s in %.1f s%n", id,
					size, p.getProperty("db"), (System.nanoTime() - start) / 1e9);
		} catch (IOException | Blastp.UsageException | RuntimeException e) {
			try {
				Files.deleteIfExists(part);
				Files.writeString(error(id), String.valueOf(e.getMessage()) + "\n", StandardCharsets.UTF_8);
			} catch (IOException e2) {
				throw new UncheckedIOException(e2);
			}
			log.println("blastp4j daemon: job " + id + " failed: " + e);
		}
		try {
			Files.deleteIfExists(fasta(id));
			Files.deleteIfExists(options(id));
		} catch (IOException e) {
			log.println("blastp4j daemon: job " + id + ": cannot delete its input: " + e.getMessage());
		}
	}

	/** The job's options, checked as blastp4j's command line checks them. */
	private Blastp.Options options(Path query, Database d, Map<String, String> params) throws Blastp.UsageException {
		List<String> args = new ArrayList<>(List.of("-query", query.toString(), "-db", d.path().toString(),
				"-num_threads", Integer.toString(threads)));
		for (Map.Entry<String, String> e : params.entrySet()) {
			if (!e.getKey().equals("db")) {
				if (!JOB_OPTIONS.contains(e.getKey())) {
					throw new Blastp.UsageException("unknown or unsupported parameter " + e.getKey());
				}
				args.add("-" + e.getKey());
				args.add(e.getValue());
			}
		}
		return Blastp.parse(args.toArray(String[]::new));
	}

	private void databases(HttpExchange x) throws IOException {
		if (!x.getRequestMethod().equals("GET")) {
			send(x, 405, "use GET\n");
			return;
		}
		StringBuilder b = new StringBuilder();
		databases.forEach((name, d) -> b.append(name).append('\t').append(d.db().size()).append('\t')
				.append(d.db().totalResidues()).append('\n'));
		send(x, 200, b.toString());
	}

	private void jobs(HttpExchange x) throws IOException {
		String path = x.getRequestURI().getPath();
		if (path.equals("/jobs") || path.equals("/jobs/")) {
			if (x.getRequestMethod().equals("POST")) {
				submit(x);
			} else {
				send(x, 405, "use POST to submit a job\n");
			}
			return;
		}
		Matcher m = JOB_PATH.matcher(path);
		if (!m.matches()) {
			send(x, 404, "no such resource\n");
			return;
		}
		if (!x.getRequestMethod().equals("GET")) {
			send(x, 405, "use GET\n");
			return;
		}
		long id;
		try {
			id = Long.parseLong(m.group(1));
		} catch (NumberFormatException e) {
			send(x, 404, "no such job\n");
			return;
		}
		State state = pending.get(id);
		if (state != null) {
			send(x, 202, state.name().toLowerCase(Locale.ROOT) + "\n");
		} else if (Files.exists(result(id))) {
			x.getResponseHeaders().set("Content-Type", "text/tab-separated-values");
			x.sendResponseHeaders(200, Files.size(result(id)));
			try (OutputStream out = x.getResponseBody()) {
				Files.copy(result(id), out);
			}
		} else if (Files.exists(error(id))) {
			send(x, 500, Files.readString(error(id), StandardCharsets.UTF_8));
		} else {
			send(x, 404, "no such job\n");
		}
	}

	/** Stores a posted job and queues it; answers with its id. */
	private void submit(HttpExchange x) throws IOException {
		Map<String, String> params = new LinkedHashMap<>();
		String query = x.getRequestURI().getRawQuery();
		if (query != null && !query.isEmpty()) {
			for (String pair : query.split("&")) {
				int eq = pair.indexOf('=');
				String name = URLDecoder.decode(eq < 0 ? pair : pair.substring(0, eq), StandardCharsets.UTF_8);
				String value = eq < 0 ? "" : URLDecoder.decode(pair.substring(eq + 1), StandardCharsets.UTF_8);
				if (params.put(name, value) != null) {
					send(x, 400, "parameter " + name + " is given twice\n");
					return;
				}
			}
		}
		Database d = databases.get(params.get("db"));
		if (d == null) {
			send(x, 400, "name a database with db=: one of " + String.join(", ", databases.keySet()) + "\n");
			return;
		}
		long id = nextId.getAndIncrement();
		Path fasta = fasta(id);
		try {
			try (InputStream in = x.getRequestBody(); OutputStream out = Files.newOutputStream(fasta)) {
				byte[] buffer = new byte[1 << 16];
				long total = 0;
				for (int n; (n = in.read(buffer)) > 0;) {
					total += n;
					if (total > maxQueryBytes) {
						throw new Blastp.UsageException("the query FASTA is larger than " + maxQueryBytes + " bytes");
					}
					out.write(buffer, 0, n);
				}
			}
			options(fasta, d, params);
			try (SequenceDatabase queries = FastaReader.read(fasta)) {
				if (queries.size() == 0) {
					throw new Blastp.UsageException("the body holds no FASTA sequences");
				}
			}
			// The options file last, under its final name only once complete: a job is in the input directory when
			// its options file is.
			Properties p = new Properties();
			params.forEach(p::setProperty);
			Path part = input.resolve(id + ".options.part");
			try (Writer w = new BufferedWriter(Files.newBufferedWriter(part, StandardCharsets.UTF_8))) {
				p.store(w, null);
			}
			Files.move(part, options(id), StandardCopyOption.ATOMIC_MOVE);
		} catch (Blastp.UsageException | IOException | RuntimeException e) {
			Files.deleteIfExists(fasta);
			send(x, 400, e.getMessage() + "\n");
			return;
		}
		pending.put(id, State.QUEUED);
		queue.add(id);
		x.getResponseHeaders().set("Location", "/jobs/" + id);
		send(x, 201, id + "\n");
	}

	private static void send(HttpExchange x, int status, String body) throws IOException {
		byte[] bytes = body.getBytes(StandardCharsets.UTF_8);
		x.getResponseHeaders().set("Content-Type", "text/plain; charset=utf-8");
		x.sendResponseHeaders(status, bytes.length);
		try (OutputStream out = x.getResponseBody()) {
			out.write(bytes);
		}
	}

	/** Stops the server and the worker, and closes the databases; queued jobs stay in the input directory. */
	@Override
	public void close() throws InterruptedException {
		server.stop(0);
		worker.interrupt();
		worker.join();
		pool.shutdownNow();
		for (Database d : databases.values()) {
			d.db().close();
		}
	}

	/** Lines written to {@code out}, each with {@code prefix}. */
	private static final class PrefixStream extends PrintStream {
		PrefixStream(PrintStream out, String prefix) {
			super(new OutputStream() {
				private boolean lineStart = true;

				@Override
				public void write(int b) {
					if (lineStart) {
						out.print(prefix);
					}
					out.write(b);
					lineStart = b == '\n';
				}
			}, true, StandardCharsets.UTF_8);
		}
	}

	/** {@code blastp4j daemon ARGS}: runs until the process is stopped; returns the exit status of a usage error. */
	static int run(String[] args, PrintStream stdout, PrintStream stderr) {
		if (args.length == 0 || List.of(args).contains("-h") || List.of(args).contains("-help")) {
			stdout.print(USAGE);
			return args.length == 0 ? 1 : 0;
		}
		Path dir = null;
		Map<String, Path> dbs = new LinkedHashMap<>();
		int port = 8080;
		String bind = "127.0.0.1";
		int threads = Runtime.getRuntime().availableProcessors();
		long maxQueryMb = 64;
		try {
			for (int i = 0; i < args.length; i += 2) {
				if (i + 1 >= args.length) {
					throw new Blastp.UsageException("option " + args[i] + " needs a value");
				}
				String v = args[i + 1];
				switch (args[i]) {
				case "-dir" -> dir = Path.of(v);
				case "-db" -> {
					int eq = v.indexOf('=');
					if (eq <= 0 || dbs.put(v.substring(0, eq), Path.of(v.substring(eq + 1))) != null) {
						throw new Blastp.UsageException("-db needs NAME=DB, each name once: '" + v + "'");
					}
				}
				case "-port" -> port = positive("-port", v, 0);
				case "-bind" -> bind = v;
				case "-num_threads" -> threads = positive("-num_threads", v, 1);
				case "-max_query_mb" -> maxQueryMb = positive("-max_query_mb", v, 1);
				default -> throw new Blastp.UsageException("unknown option " + args[i]);
				}
			}
			if (dir == null || dbs.isEmpty()) {
				throw new Blastp.UsageException("-dir and at least one -db are required");
			}
		} catch (Blastp.UsageException e) {
			stderr.println("blastp4j daemon: " + e.getMessage());
			stderr.println("Try 'blastp4j daemon -help'.");
			return 1;
		}
		try {
			Daemon daemon = new Daemon(dir, dbs, new InetSocketAddress(InetAddress.getByName(bind), port), threads,
					maxQueryMb << 20, stderr);
			Runtime.getRuntime().addShutdownHook(new Thread(() -> {
				try {
					daemon.close();
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}));
			daemon.worker.join();
			return 0;
		} catch (IOException e) {
			stderr.println("blastp4j daemon: " + e.getMessage());
			return 2;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			return 3;
		}
	}

	private static int positive(String name, String v, int min) throws Blastp.UsageException {
		try {
			int n = Integer.parseInt(v);
			if (n >= min) {
				return n;
			}
		} catch (NumberFormatException e) {
			// below
		}
		throw new Blastp.UsageException(name + " needs a whole number of at least " + min + ", not '" + v + "'");
	}
}
