package swiss.sib.blastp4j.search;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import swiss.sib.blastp4j.align.GappedAligner;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.align.HspTree;
import swiss.sib.blastp4j.align.SemiGappedLanes;
import swiss.sib.blastp4j.compo.CompoHeap;
import swiss.sib.blastp4j.compo.CompoRedo;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.seed.BatchWordFinder;
import swiss.sib.blastp4j.seed.LookupTable;
import swiss.sib.blastp4j.seed.SeedParameters;
import swiss.sib.blastp4j.seed.UngappedHsp;
import swiss.sib.blastp4j.seed.WordFinder;
import swiss.sib.blastp4j.seed.WordScan;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * A blastp search, without composition-based statistics ({@code -comp_based_stats 0}) or with mode 2 (blastp's
 * default), following BLAST+'s engine:
 * <ol>
 * <li>seeding ({@link WordFinder}) and, per subject with seeds, the preliminary gapped extensions
 * ({@code BLAST_GetGappedScore}): score only, X-drop 15 bits, from the best window of each ungapped alignment not
 * already inside a gapped one; the restricted (approximate) alignment when the best ungapped score is below 0.68 times
 * the cutoff, redone exactly when its result is inconclusive;</li>
 * <li>HSPs with common start or end points purged, E-values, the E-value cutoff; the best 550 subjects kept;</li>
 * <li>for those, the traceback ({@code Blast_TracebackFromHSPList}): each HSP not inside an already traced one is
 * aligned again from its start point with X-drop 25 bits and traceback; purging, containment, E-values again;</li>
 * <li>the best 500 subjects, ordered by E-value.</li>
 * </ol>
 * Instances are immutable; {@link #search} may be called from several threads.
 */
public final class BlastpSearch {
	/** A subject with its alignments, best first. */
	public record Hit(int subject, List<Hsp> hsps, double bestEvalue) {
	}

	/**
	 * The hits of a query, and whether the composition-based stage's heap ended filled to its cutoff
	 * ({@code BlastCompo_HeapFilledToCutoff}; always false without composition-based statistics).
	 */
	public record Result(List<Hit> hits, boolean filledToCutoff) {
	}

	/** E-value under which composition-based statistics keeps all subjects ({@code PSI_INCLUSION_ETHRESH}). */
	static final double INCLUSION_ETHRESH = 0.002;
	/** Largest E-value cutoff for which the restricted alignment may be used. */
	static final double RESTRICTED_ALIGNMENT_WORST_EVALUE = 10.0;
	/** Fraction of the cutoff below which the restricted alignment is used. */
	static final double RESTRICTED_MULT = 0.68;
	/**
	 * The "shortest subject length" BLAST uses for its cutoff scores. A BLAST database does not record it, so SeqDB
	 * reports 0 and BLAST falls back to {@code BLAST_SEQSRC_MINLENGTH}, 10, whatever the real shortest sequence is.
	 */
	public static final int MIN_SUBJECT_LENGTH = 10;
	private static final int CHUNK = 2048;
	private static final double LN2 = 0.69314718055994530941723212145818;

	private final SequenceDatabase db;
	private final Statistics statistics;
	private final int minLength;
	private final double evalue, prelimEvalue;
	private final boolean compositionBased, restrictedAlign;
	private final int hitlistSize, prelimHitlistSize, maxHsps;
	private final int gapOpen, gapExtend, xDrop, xDropFinal;

	/**
	 * A search without composition-based statistics.
	 *
	 * @param evalue the E-value cutoff (10)
	 */
	public BlastpSearch(SequenceDatabase db, Statistics statistics, double evalue) {
		this(db, statistics, evalue, false);
	}

	/**
	 * @param evalue           the E-value cutoff (10)
	 * @param compositionBased with composition-based statistics mode 2 (blastp's default)
	 */
	public BlastpSearch(SequenceDatabase db, Statistics statistics, double evalue, boolean compositionBased) {
		this(db, statistics, evalue, compositionBased, 500, 0);
	}

	/**
	 * @param evalue           the E-value cutoff ({@code -evalue}, 10)
	 * @param compositionBased with composition-based statistics mode 2 ({@code -comp_based_stats 2}, the default)
	 * @param maxTargetSeqs    number of subjects to report ({@code -max_target_seqs}, 500)
	 * @param maxHsps          HSPs to report per subject ({@code -max_hsps}); 0 for all
	 */
	public BlastpSearch(SequenceDatabase db, Statistics statistics, double evalue, boolean compositionBased,
			int maxTargetSeqs, int maxHsps) {
		if (maxTargetSeqs < 1 || maxHsps < 0) {
			throw new IllegalArgumentException("max_target_seqs must be at least 1 and max_hsps at least 0");
		}
		this.hitlistSize = maxTargetSeqs;
		this.maxHsps = maxHsps;
		this.prelimHitlistSize = prelimHitlistSize(maxTargetSeqs, compositionBased);
		this.restrictedAlign = evalue <= RESTRICTED_ALIGNMENT_WORST_EVALUE;
		this.db = db;
		this.statistics = statistics;
		this.evalue = evalue;
		this.compositionBased = compositionBased;
		// With composition-based statistics mode 2, the preliminary stage uses 5 times the E-value cutoff.
		this.prelimEvalue = compositionBased ? 5 * evalue : evalue;
		this.minLength = MIN_SUBJECT_LENGTH;
		double lambda = statistics.parameters().lambda();
		this.gapOpen = statistics.parameters().gapOpen();
		this.gapExtend = statistics.parameters().gapExtend();
		this.xDrop = (int) (15 * LN2 / lambda);
		this.xDropFinal = (int) Math.max(25 * LN2 / lambda, xDrop);
	}

	/** {@code GetPrelimHitlistSize}: how many subjects the preliminary stage keeps. */
	static int prelimHitlistSize(int hitlistSize, boolean compositionBased) {
		if (compositionBased) {
			return hitlistSize <= 500 ? 1050 : hitlistSize * 2 + 50;
		}
		return Math.min(Math.max(2 * hitlistSize, 10), hitlistSize + 50);
	}

	/** The seeding parameters this search uses for a query. */
	public SeedParameters seedParameters(byte[] query) {
		return SeedParameters.blastp(query, statistics, prelimEvalue, minLength);
	}

	/** The raw score an HSP needs to be kept ({@code cutoff_score} of BLAST's hit saving parameters). */
	public int cutoff(byte[] query) {
		return statistics.score(prelimEvalue, query.length, minLength);
	}

	/** BLAST's ordering of subjects ({@code s_EvalueCompareHSPLists}). */
	static final Comparator<Hit> BY_EVALUE = (h1, h2) -> {
		int r = compareEvalues(h1.bestEvalue, h2.bestEvalue);
		if (r == 0) {
			r = Integer.compare(h2.hsps.get(0).score, h1.hsps.get(0).score);
		}
		if (r == 0) {
			r = Integer.compare(h2.subject, h1.subject);
		}
		return r;
	};

	/** BLAST's order of the HSPs of a subject in a full hit list ({@code s_EvalueCompareHSPs}). */
	static final Comparator<Hsp> HSP_BY_EVALUE = (h1, h2) -> {
		int r = compareEvalues(h1.evalue, h2.evalue);
		return r != 0 ? r : Hsp.BY_SCORE.compare(h1, h2);
	};

	/**
	 * The final hit list: with {@code -max_hsps}, the first HSPs of each subject ({@code Blast_TrimHSPListByMaxHsps});
	 * the best {@code max_target_seqs} subjects in BLAST's order ({@code Blast_HSPResultsSortByEvalue} and pruning).
	 * Then, as BLAST's C++ layer does when it converts the results to alignments ({@code Blast_HSPListSortByEvalue} in
	 * {@code blast_seqalign.cpp}), the HSPs of each subject are sorted by E-value and score; with E-values below
	 * 1e-180 counted as equal, HSPs of equal (final) score are then ordered by position.
	 */
	List<Hit> finalHitList(List<Hit> hits) {
		if (maxHsps > 0) {
			for (Hit h : hits) {
				if (h.hsps().size() > maxHsps) {
					h.hsps().subList(maxHsps, h.hsps().size()).clear();
				}
			}
		}
		hits.sort(BY_EVALUE);
		List<Hit> kept = new ArrayList<>(hits.subList(0, Math.min(hitlistSize, hits.size())));
		for (Hit h : kept) {
			h.hsps().sort(HSP_BY_EVALUE);
		}
		return kept;
	}

	/** {@code s_EvalueComp}: E-values below 1e-180 are equal. */
	static int compareEvalues(double e1, double e2) {
		if (e1 < 1.0e-180 && e2 < 1.0e-180) {
			return 0;
		}
		return Double.compare(e1, e2) < 0 ? -1 : e1 > e2 ? 1 : 0;
	}

	/**
	 * Searches a query (NCBIstdaa codes) with {@code tasks} tasks of the pool, as BLAST does a query that is alone in
	 * its batch.
	 */
	public List<Hit> search(byte[] query, ExecutorService pool, int tasks) throws InterruptedException {
		return search(query, pool, tasks, true).hits();
	}

	/**
	 * Searches a query (NCBIstdaa codes) with {@code tasks} tasks of the pool.
	 *
	 * @param earlyTermination whether the composition-based stage may skip subjects once its heap is filled to the
	 *                         cutoff; BLAST does so only when the heaps of all queries in the batch are (see
	 *                         {@link QueryBatches})
	 */
	public Result search(byte[] query, ExecutorService pool, int tasks, boolean earlyTermination)
			throws InterruptedException {
		SeedParameters seeds = seedParameters(query);
		LookupTable table = LookupTable.of(query, seeds.threshold());
		int cutoff = cutoff(query);

		// Preliminary stage.
		AtomicInteger next = new AtomicInteger();
		List<Callable<List<Hit>>> prelimTasks = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			prelimTasks.add(() -> prelim(query, WordFinder.of(query, table, seeds), cutoff, next));
		}
		List<Hit> prelim = new ArrayList<>();
		for (List<Hit> part : runAll(pool, prelimTasks)) {
			prelim.addAll(part);
		}
		prelim.sort(BY_EVALUE);
		return finish(query, prelim.subList(0, Math.min(prelimHitlistSize, prelim.size())), pool, tasks,
				earlyTermination);
	}

	/**
	 * The preliminary stage for a batch of queries, as BLAST runs it for a batch: each subject is scanned once for all
	 * of them ({@link BatchWordFinder}). Gives each query's preliminary hits, best first and cut to the preliminary
	 * hit list size, the same as {@link #search} keeps for the query alone; {@link #finish} completes each query.
	 */
	public List<List<Hit>> prelim(List<byte[]> queries, ExecutorService pool, int tasks) throws InterruptedException {
		final int n = queries.size();
		final byte[][] qs = queries.toArray(byte[][]::new);
		final SeedParameters[] seeds = new SeedParameters[n];
		final LookupTable[] tables = new LookupTable[n];
		final int[] cutoffs = new int[n];
		for (int i = 0; i < n; i++) {
			seeds[i] = seedParameters(qs[i]);
			tables[i] = LookupTable.of(qs[i], seeds[i].threshold());
			cutoffs[i] = cutoff(qs[i]);
		}
		final BatchWordFinder.Table table = BatchWordFinder.Table.of(tables);
		AtomicInteger next = new AtomicInteger();
		List<Callable<List<List<Hit>>>> prelimTasks = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			prelimTasks.add(() -> prelim(qs, new BatchWordFinder(qs, tables, seeds, table), cutoffs, next));
		}
		List<List<Hit>> found = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			found.add(new ArrayList<>());
		}
		for (List<List<Hit>> part : runAll(pool, prelimTasks)) {
			for (int i = 0; i < n; i++) {
				found.get(i).addAll(part.get(i));
			}
		}
		List<List<Hit>> kept = new ArrayList<>();
		for (List<Hit> prelim : found) {
			prelim.sort(BY_EVALUE);
			kept.add(prelim.subList(0, Math.min(prelimHitlistSize, prelim.size())));
		}
		return kept;
	}

	/**
	 * The stages after the preliminary one for a query: composition-based statistics or the traceback, as
	 * {@link #search} does with its preliminary hits.
	 *
	 * @param kept             the query's preliminary hits, best first, as {@link #prelim(List, ExecutorService, int)}
	 *                         gives them
	 * @param earlyTermination as for {@link #search(byte[], ExecutorService, int, boolean)}
	 */
	public Result finish(byte[] query, List<Hit> kept, ExecutorService pool, int tasks, boolean earlyTermination)
			throws InterruptedException {
		if (compositionBased) {
			return compositionBasedStage(query, kept, pool, tasks, earlyTermination);
		}

		// Traceback stage.
		AtomicInteger nextHit = new AtomicInteger();
		List<Callable<List<Hit>>> tracebackTasks = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			tracebackTasks.add(() -> {
				GappedAligner aligner = new GappedAligner(gapOpen, gapExtend);
				List<Hit> out = new ArrayList<>();
				byte[] subject = null;
				for (int i; (i = nextHit.getAndIncrement()) < kept.size();) {
					subject = db.withSentinels(kept.get(i).subject(), subject);
					Hit h = traceback(query, aligner, kept.get(i), subject);
					if (h != null) {
						out.add(h);
					}
				}
				return out;
			});
		}
		List<Hit> hits = new ArrayList<>();
		for (List<Hit> part : runAll(pool, tracebackTasks)) {
			hits.addAll(part);
		}
		return new Result(finalHitList(hits), false);
	}

	/**
	 * The traceback stage with composition-based statistics ({@code Blast_RedoAlignmentCore_MT}). BLAST reads the
	 * preliminary subjects best first (the HSP stream sorts on read: {@code Blast_HSPResultsReverseSort}) and skips a
	 * subject when its result heap is already full of subjects below the inclusion E-value and the subject's
	 * preliminary E-value is more than 5 times that. The redoing is done in parallel here; the skipping and the heap
	 * are then replayed in BLAST's order, which gives the same result because redoing one subject does not depend on
	 * another.
	 * <p>
	 * This is BLAST's result with one thread. With more threads BLAST gives each thread its own heap and skipping over
	 * a static share of the subjects (of all queries in the batch), so its output can depend on the number of threads
	 * when a query has more than {@code max_target_seqs} significant subjects.
	 *
	 * @param earlyTermination whether subjects may be skipped at all (see {@link #search(byte[], ExecutorService, int,
	 *                         boolean)})
	 */
	private Result compositionBasedStage(byte[] query, List<Hit> prelim, ExecutorService pool, int tasks,
			boolean earlyTermination) throws InterruptedException {
		CompoRedo.Query q = new CompoRedo.Query(query);
		List<Hit> ordered = new ArrayList<>(prelim);
		ordered.sort(BY_EVALUE);
		CompoRedo.Result[] results = new CompoRedo.Result[ordered.size()];
		AtomicInteger next = new AtomicInteger();
		List<Callable<Void>> work = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			work.add(() -> {
				CompoRedo redo = new CompoRedo(q, statistics);
				byte[] subject = null;
				for (int i; (i = next.getAndIncrement()) < ordered.size();) {
					Hit h = ordered.get(i);
					subject = db.withSentinels(h.subject(), subject);
					results[i] = redo.redo(h.hsps(), subject, 1, db.length(h.subject()), evalue, null);
				}
				return null;
			});
		}
		runAll(pool, work);

		CompoHeap<Hit> heap = new CompoHeap<>(hitlistSize, INCLUSION_ETHRESH);
		// BLAST's thread-local matrix: unscaled BLOSUM62 until the first successful adjustment.
		int[][] current = new int[swiss.sib.blastp4j.stats.Blosum62.SIZE][swiss.sib.blastp4j.stats.Blosum62.SIZE];
		for (int a = 0; a < current.length; a++) {
			for (int b = 0; b < current.length; b++) {
				current[a][b] = swiss.sib.blastp4j.stats.Blosum62.score(a, b);
			}
		}
		CompoRedo sequential = null;
		for (int i = 0; i < ordered.size(); i++) {
			Hit h = ordered.get(i);
			if (earlyTermination && heap.earlyTermination(h.bestEvalue())) {
				continue;
			}
			CompoRedo.Result r = results[i];
			if (r.usedInheritedMatrix()) {
				// Rare: redo this subject with the matrix BLAST would have at this point.
				if (sequential == null) {
					sequential = new CompoRedo(q, statistics);
				}
				r = sequential.redo(h.hsps(), db.withSentinels(h.subject(), null), 1, db.length(h.subject()), evalue,
						current);
			}
			if (r.finalMatrix() != null) {
				current = r.finalMatrix();
			}
			if (r.hsps().isEmpty() || !(r.bestEvalue() <= evalue)) {
				continue;
			}
			if (heap.wouldInsert(r.bestEvalue(), r.bestScaledScore(), h.subject())) {
				heap.insert(new Hit(h.subject(), r.hsps(), r.bestEvalue()), r.bestEvalue(), r.bestScaledScore(),
						h.subject());
			}
		}
		boolean filled = heap.filledToCutoff();
		List<Hit> hits = new ArrayList<>();
		for (Hit h; (h = heap.pop()) != null;) {
			hits.add(h);
		}
		return new Result(finalHitList(hits), filled);
	}

	private static <T> List<T> runAll(ExecutorService pool, List<Callable<T>> tasks) throws InterruptedException {
		List<Future<T>> futures = new ArrayList<>();
		for (Callable<T> t : tasks) {
			futures.add(pool.submit(t));
		}
		List<T> out = new ArrayList<>();
		try {
			for (Future<T> f : futures) {
				out.add(f.get());
			}
		} catch (ExecutionException e) {
			throw new IllegalStateException(e.getCause());
		}
		return out;
	}

	/** BLAST's order of the ungapped alignments ({@code score_compare_match}). */
	static final Comparator<UngappedHsp> SEED_ORDER = (h1, h2) -> {
		int r = Integer.compare(h2.score(), h1.score());
		if (r == 0) {
			r = Integer.compare(h1.subjectStart(), h2.subjectStart());
		}
		if (r == 0) {
			r = Integer.compare(h2.length(), h1.length());
		}
		if (r == 0) {
			r = Integer.compare(h1.queryStart(), h2.queryStart());
		}
		return r;
	};

	List<Hit> prelim(byte[] query, WordScan finder, int cutoff, AtomicInteger next) {
		final GappedAligner aligner = new GappedAligner(gapOpen, gapExtend);
		final List<Hit> found = new ArrayList<>();
		final List<UngappedHsp> seeds = new ArrayList<>();
		final LaneStage lanes = SemiGappedLanes.AVAILABLE ? new LaneStage(aligner, query.length, List.of(found)) : null;
		byte[] subject = null;
		for (int chunk; (chunk = next.getAndIncrement()) * CHUNK < db.size();) {
			int end = Math.min(db.size(), (chunk + 1) * CHUNK);
			for (int i = chunk * CHUNK; i < end; i++) {
				subject = db.withSentinels(i, subject);
				int length = db.length(i);
				if (finder.vectorScan()) {
					finder.findBlocks(subject, 1, length, seeds);
				} else {
					finder.find(subject, 1, length, seeds);
				}
				if (seeds.isEmpty()) {
					continue;
				}
				if (lanes != null) {
					lanes.add(0, query, cutoff, i, Arrays.copyOf(subject, length + 2), length, seeds);
					seeds.clear();
					continue;
				}
				Hit hit = prelimHit(query, aligner, i, subject, length, seeds, cutoff);
				seeds.clear();
				if (hit != null) {
					found.add(hit);
				}
			}
		}
		if (lanes != null) {
			lanes.finish();
		}
		return found;
	}

	/** The preliminary stage of a batch's task: {@link #prelim(byte[], WordScan, int, AtomicInteger)} per query. */
	@SuppressWarnings("unchecked")
	private List<List<Hit>> prelim(byte[][] queries, BatchWordFinder finder, int[] cutoffs, AtomicInteger next) {
		final GappedAligner aligner = new GappedAligner(gapOpen, gapExtend);
		final int n = queries.length;
		final List<List<Hit>> found = new ArrayList<>();
		final List<UngappedHsp>[] seeds = new List[n];
		int longest = 0;
		for (int q = 0; q < n; q++) {
			found.add(new ArrayList<>());
			seeds[q] = new ArrayList<>();
			longest = Math.max(longest, queries[q].length);
		}
		final LaneStage lanes = SemiGappedLanes.AVAILABLE ? new LaneStage(aligner, longest, found) : null;
		byte[] subject = null;
		for (int chunk; (chunk = next.getAndIncrement()) * CHUNK < db.size();) {
			int end = Math.min(db.size(), (chunk + 1) * CHUNK);
			for (int i = chunk * CHUNK; i < end; i++) {
				subject = db.withSentinels(i, subject);
				int length = db.length(i);
				finder.find(subject, 1, length, seeds);
				// The lanes' jobs keep the subject: one copy for all queries.
				byte[] residues = null;
				for (int q = 0; q < n; q++) {
					if (seeds[q].isEmpty()) {
						continue;
					}
					if (lanes != null) {
						if (residues == null) {
							residues = Arrays.copyOf(subject, length + 2);
						}
						lanes.add(q, queries[q], cutoffs[q], i, residues, length, seeds[q]);
						seeds[q].clear();
						continue;
					}
					Hit hit = prelimHit(queries[q], aligner, i, subject, length, seeds[q], cutoffs[q]);
					seeds[q].clear();
					if (hit != null) {
						found.get(q).add(hit);
					}
				}
			}
		}
		if (lanes != null) {
			lanes.finish();
		}
		return found;
	}

	/**
	 * The gapped scores of the preliminary stage ({@link #gappedScores}) with the score-only extensions in vector lanes
	 * ({@link SemiGappedLanes}). The extensions of one subject depend on each other (a seed within an HSP found
	 * before is skipped), so the lanes take those of many subjects: each subject is a job that runs its seeds as
	 * {@link #gappedScores} does until a seed needs its two extensions, submits them, and continues when both are
	 * done. The restricted alignment stays scalar. Each query's hits are as {@link #prelimHit} gives them, in subject
	 * order.
	 */
	private final class LaneStage implements SemiGappedLanes.Done {
		/** One subject's seeds for one query, as {@link #gappedScores} runs them. */
		private final class Job {
			final int q, subject, length, cutoff, restrictedCutoff;
			final byte[] query, residues;
			final List<UngappedHsp> seeds;
			final List<Hsp> hsps = new ArrayList<>();
			final HspTree tree = new HspTree();
			boolean restricted;
			int index, slot;
			// The seed being extended: its start point, extensions pending, and their results.
			int qOff, sOff, pending, score, queryStart, queryStop, subjectStart, subjectStop;

			Job(int q, byte[] query, int cutoff, int subject, byte[] residues, int length, List<UngappedHsp> seeds) {
				this.q = q;
				this.query = query;
				this.cutoff = cutoff;
				this.subject = subject;
				this.residues = residues;
				this.length = length;
				this.seeds = seeds;
				seeds.sort(SEED_ORDER);
				restrictedCutoff = (int) (RESTRICTED_MULT * cutoff);
				restricted = restrictedAlign && seeds.get(0).score() < restrictedCutoff;
			}
		}

		private final GappedAligner aligner;
		private final SemiGappedLanes lanes;
		private final List<List<Hit>> found;
		private Job[] jobs = new Job[64];
		private int[] free = new int[64];
		private int freeCount;
		private int used;

		LaneStage(GappedAligner aligner, int longestQuery, List<List<Hit>> found) {
			this.aligner = aligner;
			this.lanes = new SemiGappedLanes(gapOpen, gapExtend, longestQuery);
			this.found = found;
		}

		/**
		 * Adds the seeds of query {@code q} in subject {@code i} ({@code residues}, which the job keeps, as
		 * {@link SequenceDatabase#withSentinels} gives it), then runs the lanes while enough is queued.
		 */
		void add(int q, byte[] query, int cutoff, int i, byte[] residues, int length, List<UngappedHsp> seeds) {
			final Job job = new Job(q, query, cutoff, i, residues, length, new ArrayList<>(seeds));
			job.slot = slot(job);
			run(job);
			while (lanes.queued() >= 2 * SemiGappedLanes.L) {
				lanes.step(this);
			}
		}

		/** Runs the lanes until every job is done, then puts each query's hits in subject order. */
		void finish() {
			while (lanes.step(this)) {
			}
			for (List<Hit> hits : found) {
				hits.sort(Comparator.comparingInt(Hit::subject));
			}
		}

		private int slot(Job job) {
			if (freeCount > 0) {
				final int s = free[--freeCount];
				jobs[s] = job;
				return s;
			}
			if (used == jobs.length) {
				jobs = Arrays.copyOf(jobs, 2 * used);
				free = Arrays.copyOf(free, 2 * used);
			}
			jobs[used] = job;
			return used++;
		}

		/** {@link #gappedScores}'s loop from the job's seed on, until a seed's extensions are submitted. */
		private void run(Job j) {
			final int base = 1;
			for (; j.index < j.seeds.size(); j.index++) {
				final UngappedHsp seed = j.seeds.get(j.index);
				final int qStart = seed.queryStart(), sStart = seed.subjectStart();
				final int qEnd = qStart + seed.length(), sEnd = sStart + seed.length();
				if (j.tree.contains(seed.score(), qStart, qEnd, sStart, sEnd)) {
					continue;
				}
				final int qOff = aligner.startForGapped(j.query, j.residues, base, qStart, seed.length(), sStart,
						seed.length());
				final int sOff = seed.subjectSeed() + qOff - seed.querySeed();
				if (j.restricted) {
					aligner.extendScoreOnly(j.query, j.residues, base, j.length, qOff, sOff, true, xDrop);
					if (aligner.score < j.cutoff && aligner.score >= j.restrictedCutoff) {
						// Inconclusive: start again with exact alignments for this subject.
						j.restricted = false;
						j.hsps.clear();
						j.tree.clear();
						j.index = -1;
						continue;
					}
					if (aligner.score >= j.cutoff) {
						final Hsp h = new Hsp(aligner.score, aligner.queryStart, aligner.queryStop,
								aligner.subjectStart, aligner.subjectStop, qOff, sOff);
						j.hsps.add(h);
						j.tree.add(h);
					}
					continue;
				}
				// GappedAligner.extendScoreOnly's two extensions: left including the start point, then right.
				final int qLength = qOff + 1, sLength = sOff + 1;
				j.qOff = qOff;
				j.sOff = sOff;
				j.score = 0;
				j.queryStop = qOff;
				j.subjectStop = sOff;
				j.pending = 1;
				lanes.submit(j.query, 0, qLength, j.residues, base, sLength, true, xDrop, 2 * j.slot);
				if (qLength < j.query.length && sLength < j.length) {
					j.pending = 2;
					lanes.submit(j.query, qOff, j.query.length - qLength, j.residues, base + sOff,
							j.length - sLength, false, xDrop, 2 * j.slot + 1);
				}
				return;
			}
			// All seeds done.
			jobs[j.slot] = null;
			free[freeCount++] = j.slot;
			final Hit hit = prelimHit(j.query, j.subject, j.length, j.hsps);
			if (hit != null) {
				found.get(j.q).add(hit);
			}
		}

		@Override
		public void done(int id, int score, int aOffset, int bOffset) {
			final Job j = jobs[id >> 1];
			if ((id & 1) == 0) {
				j.queryStart = j.qOff + 1 - aOffset;
				j.subjectStart = j.sOff + 1 - bOffset;
			} else {
				j.queryStop = aOffset + j.qOff + 1;
				j.subjectStop = bOffset + j.sOff + 1;
			}
			j.score += score;
			if (--j.pending > 0) {
				return;
			}
			if (j.score >= j.cutoff) {
				final Hsp h = new Hsp(j.score, j.queryStart, j.queryStop, j.subjectStart, j.subjectStop, j.qOff,
						j.sOff);
				j.hsps.add(h);
				j.tree.add(h);
			}
			j.index++;
			run(j);
		}
	}

	/**
	 * The gapped scores of the query's seeds in subject {@code i} and the HSPs the preliminary stage keeps of them, as
	 * a hit; null when it keeps none.
	 */
	private Hit prelimHit(byte[] query, GappedAligner aligner, int i, byte[] subject, int length,
			List<UngappedHsp> seeds, int cutoff) {
		return prelimHit(query, i, length, gappedScores(query, aligner, subject, length, seeds, cutoff));
	}

	/** The HSPs the preliminary stage keeps of the gapped scores {@code hsps} of subject {@code i}, as a hit. */
	private Hit prelimHit(byte[] query, int i, int length, List<Hsp> hsps) {
		if (hsps.isEmpty()) {
			return null;
		}
		purgeCommonEndpoints(hsps);
		hsps.sort(Hsp.BY_SCORE);
		double best = Double.MAX_VALUE;
		List<Hsp> kept = new ArrayList<>();
		for (Hsp h : hsps) {
			h.evalue = statistics.evalue(h.score, query.length, length);
			if (h.evalue <= prelimEvalue) {
				kept.add(h);
				best = Math.min(best, h.evalue);
			}
		}
		return kept.isEmpty() ? null : new Hit(i, kept, best);
	}

	/**
	 * {@code BLAST_GetGappedScore} for one subject, given as {@link SequenceDatabase#withSentinels} gives it: residues
	 * at {@code [1, 1 + length)}.
	 */
	List<Hsp> gappedScores(byte[] query, GappedAligner aligner, byte[] residues, int length, List<UngappedHsp> seeds,
			int cutoff) {
		final int base = 1;
		seeds.sort(SEED_ORDER);
		final int restrictedCutoff = (int) (RESTRICTED_MULT * cutoff);
		boolean restricted = restrictedAlign && seeds.get(0).score() < restrictedCutoff;
		List<Hsp> hsps = new ArrayList<>();
		HspTree tree = new HspTree();
		for (int index = 0; index < seeds.size(); index++) {
			UngappedHsp seed = seeds.get(index);
			int qStart = seed.queryStart(), sStart = seed.subjectStart();
			int qEnd = qStart + seed.length(), sEnd = sStart + seed.length();
			if (tree.contains(seed.score(), qStart, qEnd, sStart, sEnd)) {
				continue;
			}
			int qOff = aligner.startForGapped(query, residues, base, qStart, seed.length(), sStart, seed.length());
			int sOff = seed.subjectSeed() + qOff - seed.querySeed();
			aligner.extendScoreOnly(query, residues, base, length, qOff, sOff, restricted, xDrop);
			if (restricted && aligner.score < cutoff && aligner.score >= restrictedCutoff) {
				// Inconclusive: start again with exact alignments for this subject.
				restricted = false;
				hsps.clear();
				tree.clear();
				index = -1;
				continue;
			}
			if (aligner.score >= cutoff) {
				Hsp h = new Hsp(aligner.score, aligner.queryStart, aligner.queryStop, aligner.subjectStart,
						aligner.subjectStop, qOff, sOff);
				hsps.add(h);
				tree.add(h);
			}
		}
		return hsps;
	}

	/**
	 * {@code Blast_TracebackFromHSPList} for one subject, given as {@link SequenceDatabase#withSentinels} gives it;
	 * null when no HSP is left.
	 */
	Hit traceback(byte[] query, GappedAligner aligner, Hit prelim, byte[] residues) {
		final int base = 1, length = db.length(prelim.subject());
		List<Hsp> done = new ArrayList<>();
		HspTree tree = new HspTree();
		for (Hsp p : prelim.hsps()) {
			if (tree.contains(p.score, p.queryStart, p.queryEnd, p.subjectStart, p.subjectEnd)) {
				continue;
			}
			int qStart = p.queryGappedStart, sStart = p.subjectGappedStart;
			if (qStart == 0 && sStart == 0) {
				int[] start = aligner.offsetsForGapped(query, residues, base, p.queryStart, p.queryEnd, p.subjectStart,
						p.subjectEnd);
				if (start == null) {
					continue;
				}
				qStart = start[0];
				sStart = start[1];
			}
			aligner.extendWithTraceback(query, residues, base, length, qStart, sStart, xDropFinal);
			Hsp h = new Hsp(aligner.score, aligner.queryStart, aligner.queryStop, aligner.subjectStart,
					aligner.subjectStop, qStart, sStart);
			h.editScript = aligner.editScript;
			done.add(h);
			tree.add(h);
		}
		purgeCommonEndpoints(done);
		done.sort(Hsp.BY_SCORE);
		tree.clear();
		List<Hsp> kept = new ArrayList<>();
		for (Hsp h : done) {
			if (!tree.contains(h.score, h.queryStart, h.queryEnd, h.subjectStart, h.subjectEnd)) {
				kept.add(h);
				tree.add(h);
			}
		}
		List<Hsp> out = new ArrayList<>();
		double best = Double.MAX_VALUE;
		for (Hsp h : kept) {
			h.evalue = statistics.evalue(h.score, query.length, length);
			if (h.evalue <= evalue) {
				out.add(h);
				best = Math.min(best, h.evalue);
			}
		}
		return out.isEmpty() ? null : new Hit(prelim.subject(), out, best);
	}

	/**
	 * {@code Blast_HSPListPurgeHSPsWithCommonEndpoints} with purging (as for blastp): of the HSPs with the same start
	 * (query and subject), only the first in BLAST's start order is kept, then likewise for the same end.
	 */
	static void purgeCommonEndpoints(List<Hsp> hsps) {
		if (hsps.size() < 2) {
			return;
		}
		hsps.sort(BlastpSearch::compareHsp);
		for (int i = 0; i < hsps.size(); i++) {
			Hsp a = hsps.get(i);
			while (i + 1 < hsps.size() && hsps.get(i + 1).queryStart == a.queryStart
					&& hsps.get(i + 1).subjectStart == a.subjectStart) {
				hsps.remove(i + 1);
			}
		}
		hsps.sort((h1, h2) -> {
			int r = Integer.compare(h1.queryEnd, h2.queryEnd);
			if (r == 0) {
				r = Integer.compare(h1.subjectEnd, h2.subjectEnd);
			}
			if (r == 0) {
				r = Integer.compare(h2.score, h1.score);
			}
			if (r == 0) {
				r = Integer.compare(h2.queryStart, h1.queryStart);
			}
			if (r == 0) {
				r = Integer.compare(h2.subjectStart, h1.subjectStart);
			}
			return r;
		});
		for (int i = 0; i < hsps.size(); i++) {
			Hsp a = hsps.get(i);
			while (i + 1 < hsps.size() && hsps.get(i + 1).queryEnd == a.queryEnd
					&& hsps.get(i + 1).subjectEnd == a.subjectEnd) {
				hsps.remove(i + 1);
			}
		}
	}

	private static int compareHsp(Hsp h1, Hsp h2) {
		int r = Integer.compare(h1.queryStart, h2.queryStart);
		if (r == 0) {
			r = Integer.compare(h1.subjectStart, h2.subjectStart);
		}
		if (r == 0) {
			r = Integer.compare(h2.score, h1.score);
		}
		if (r == 0) {
			r = Integer.compare(h2.queryEnd, h1.queryEnd);
		}
		if (r == 0) {
			r = Integer.compare(h2.subjectEnd, h1.subjectEnd);
		}
		return r;
	}
}
