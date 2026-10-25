package swiss.sib.blastp4j.search;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.seed.LookupTable;
import swiss.sib.blastp4j.seed.SeedParameters;
import swiss.sib.blastp4j.seed.UngappedHsp;
import swiss.sib.blastp4j.seed.WordFinder;

/**
 * The seeding stage of blastp over a database: the lookup table scan, two-hit trigger and ungapped extension for every
 * subject. The subjects with at least one ungapped alignment at the cutoff are those BLAST extends with gaps.
 * <p>
 * Instances are immutable; {@link #search} may be called from several threads.
 */
public final class SeedSearch {
	/** A subject with its ungapped alignments, ordered by score (descending) as BLAST sorts them. */
	public record Seeds(int subject, List<UngappedHsp> hsps) {
	}

	/** Subjects per work unit. */
	private static final int CHUNK = 2048;

	private final SequenceDatabase db;

	public SeedSearch(SequenceDatabase db) {
		this.db = db;
	}

	/** Seeds of a query (NCBIstdaa codes) in all subjects, ordered by subject, using {@code tasks} tasks of the pool. */
	public List<Seeds> search(byte[] query, SeedParameters parameters, ExecutorService pool, int tasks)
			throws InterruptedException {
		LookupTable table = LookupTable.of(query, parameters.threshold());
		AtomicInteger next = new AtomicInteger();
		List<Future<List<Seeds>>> parts = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			parts.add(pool.submit(() -> scan(new WordFinder(query, table, parameters), next)));
		}
		List<Seeds> all = new ArrayList<>();
		try {
			for (Future<List<Seeds>> part : parts) {
				all.addAll(part.get());
			}
		} catch (ExecutionException e) {
			throw new IllegalStateException(e.getCause());
		}
		all.sort(Comparator.comparingInt(Seeds::subject));
		return all;
	}

	private List<Seeds> scan(WordFinder finder, AtomicInteger next) {
		final byte[] residues = db.residues();
		final List<Seeds> found = new ArrayList<>();
		List<UngappedHsp> hsps = new ArrayList<>();
		for (int chunk; (chunk = next.getAndIncrement()) * CHUNK < db.size();) {
			int end = Math.min(db.size(), (chunk + 1) * CHUNK);
			for (int i = chunk * CHUNK; i < end; i++) {
				finder.find(residues, db.start(i), db.length(i), hsps);
				if (!hsps.isEmpty()) {
					hsps.sort(Comparator.comparingInt(UngappedHsp::score).reversed());
					found.add(new Seeds(i, List.copyOf(hsps)));
					hsps.clear();
				}
			}
		}
		return found;
	}
}
