package swiss.sib.blastp4j.compo;

import swiss.sib.blastp4j.align.EditScript;
import swiss.sib.blastp4j.align.GappedAligner;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.align.SemiGappedLanes;
import swiss.sib.blastp4j.align.TracebackKernel;
import swiss.sib.blastp4j.align.TracebackLanes;
import swiss.sib.blastp4j.align.TracebackPairs;

/**
 * The composition-based stage of one task with its tracebacks as lanes ({@link TracebackLanes}), or where those do
 * not pay, two at a time ({@link TracebackPairs}): up to
 * {@link #IN_FLIGHT} subjects at a time, each run in {@link CompoRedo.Subject}'s steps; the two sides of each
 * traceback (alignEx from the start point to the left and to the right) go to the lanes, and their results are joined
 * as {@link GappedAligner#extendWithTraceback} joins them. Each subject's result is {@link CompoRedo#redo}'s. A
 * subject whose matrix adjustment failed uses the aligner's matrix, as in {@link CompoRedo#redo}: its tracebacks stay
 * scalar (and its result is marked stale, for the caller to redo in order).
 */
public final class TracebackStage {
	/**
	 * The tracebacks two at a time ({@link TracebackPairs}) instead of in lanes: where the lanes do not pay (the
	 * vectors have fewer than 8 ints, as AArch64's), or with {@code -Dblastp4j.tracebackPairs=true}.
	 */
	static final boolean PAIRS = !SemiGappedLanes.AVAILABLE || Boolean.getBoolean("blastp4j.tracebackPairs");
	/** Whether the stage runs: unless {@code -Dblastp4j.tracebackLanes=false} (then the scalar redo runs). */
	public static final boolean AVAILABLE = Boolean.parseBoolean(System.getProperty("blastp4j.tracebackLanes", "true"));
	/** Subjects in flight: two sides each, so that the kernel's queue stays fuller than its lanes. */
	static final int IN_FLIGHT = 2 * (PAIRS ? 2 : TracebackLanes.L);

	/**
	 * Each thread's lanes, kept from query to query: their buffers (the traceback ring above all) would otherwise be
	 * garbage after every query, which made G1 grow the heap to twice its size. Lanes sized for a longer query do for
	 * a shorter one.
	 */
	private static final ThreadLocal<TracebackKernel> LANES = new ThreadLocal<>();

	/** A subject to redo, and its index for the caller. */
	public record Job(int index, CompoRedo.Subject subject) {
	}

	/** The subjects to redo, one at a time; null when there are no more. */
	public interface Source {
		Job next();
	}

	/** Receives each subject's result. */
	public interface Sink {
		void done(int index, CompoRedo.Result result);
	}

	private final CompoRedo redo;
	private final GappedAligner aligner;
	private final TracebackKernel lanes;
	private final byte[] q;
	private final int xDrop;
	// Per slot: the job, the matrix it was last submitted with (as the lanes take it), and the two sides' results.
	private final Job[] job = new Job[IN_FLIGHT];
	private final CompoRedo.Subject[] subject = new CompoRedo.Subject[IN_FLIGHT];
	private final Object[] matrixOf = new Object[IN_FLIGHT];
	private final int[][] lanesMatrix = new int[IN_FLIGHT][1024];
	private final long[] version = new long[IN_FLIGHT];
	/** The versions of the matrices given to the lanes: each new one a new version, above BLOSUM62's 0. */
	private static final java.util.concurrent.atomic.AtomicLong VERSIONS = new java.util.concurrent.atomic.AtomicLong();
	private final int[] waiting = new int[IN_FLIGHT];
	private final int[][] side = new int[2 * IN_FLIGHT][3];
	private final EditScript[] ops = new EditScript[2 * IN_FLIGHT];
	private final EditScript none = new EditScript();
	private int inFlight;
	private Sink sink;

	public TracebackStage(CompoRedo redo) {
		this.redo = redo;
		aligner = redo.aligner();
		q = redo.queryData();
		xDrop = redo.xDrop();
		TracebackKernel l = LANES.get();
		if (l == null || !l.fits(aligner.gapOpen(), aligner.gapExtend(), q.length)) {
			final int maxA = Math.max(q.length, l == null ? 0 : l.maxA());
			l = PAIRS ? new TracebackPairs(aligner.gapOpen(), aligner.gapExtend(), maxA)
					: new TracebackLanes(aligner.gapOpen(), aligner.gapExtend(), maxA);
			LANES.set(l);
		}
		lanes = l;
	}

	/** Redoes the subjects of {@code source}, giving each result to {@code sink}. */
	public void run(Source source, Sink sink) {
		this.sink = sink;
		boolean more = true;
		while (true) {
			for (int slot = 0; more && inFlight < IN_FLIGHT && slot < IN_FLIGHT; slot++) {
				if (job[slot] == null) {
					final Job j = source.next();
					if (j == null) {
						more = false;
					} else {
						job[slot] = j;
						subject[slot] = j.subject();
						inFlight++;
						drive(slot);
					}
				}
			}
			if (inFlight == 0) {
				if (!more) {
					return;
				}
				continue;
			}
			lanes.step(this::sideDone);
		}
	}

	/** Runs the slot's subject to its next traceback, submitted to the lanes, or to its result. */
	private void drive(int slot) {
		final CompoRedo.Subject s = subject[slot];
		while (s.advance()) {
			final Hsp p = s.traceback();
			final int length = s.subject().length - 2;
			if (s.matrix() == null) {
				// As redo: the aligner's matrix, whatever it is.
				aligner.extendWithTraceback(q, s.subject(), 1, length, p.queryGappedStart, p.subjectGappedStart,
						xDrop);
				s.tracebackDone(aligner);
				continue;
			}
			if (matrixOf[slot] != s.matrix()) {
				// The slot's array refilled under a new version: none of the slot's extensions is queued now.
				final int[] m32 = lanesMatrix[slot];
				version[slot] = VERSIONS.incrementAndGet();
				for (int a = 0; a < 28; a++) {
					for (int b = 0; b < 28; b++) {
						m32[a << 5 | b] = s.matrix().score(a, b);
					}
				}
				matrixOf[slot] = s.matrix();
			}
			final int qs = p.queryGappedStart, ss = p.subjectGappedStart;
			lanes.submit(q, 0, qs + 1, s.subject(), 1, ss + 1, true, xDrop, 2 * slot, lanesMatrix[slot],
					version[slot]);
			waiting[slot] = 1;
			final int m = q.length - qs - 1, n = length - ss - 1;
			if (qs < q.length && ss < length && m > 0 && n > 0) {
				lanes.submit(q, qs, m, s.subject(), 1 + ss, n, false, xDrop, 2 * slot + 1, lanesMatrix[slot],
						version[slot]);
				waiting[slot]++;
			} else {
				// alignEx of an empty side: 0, at 0, 0.
				side[2 * slot + 1][0] = 0;
				side[2 * slot + 1][1] = 0;
				side[2 * slot + 1][2] = 0;
				ops[2 * slot + 1] = none;
			}
			return;
		}
		sink.done(job[slot].index(), s.result());
		job[slot] = null;
		subject[slot] = null;
		matrixOf[slot] = null;
		inFlight--;
	}

	private void sideDone(int id, int score, int aOffset, int bOffset, EditScript o) {
		side[id][0] = score;
		side[id][1] = aOffset;
		side[id][2] = bOffset;
		ops[id] = o.copy();
		final int slot = id >> 1;
		if (--waiting[slot] > 0) {
			return;
		}
		final CompoRedo.Subject s = subject[slot];
		final Hsp p = s.traceback();
		final int[] l = side[2 * slot], r = side[2 * slot + 1];
		aligner.joinTraceback(q.length, s.subject().length - 2, p.queryGappedStart, p.subjectGappedStart, l[0], l[1],
				l[2], ops[2 * slot], r[0], r[1], r[2], ops[2 * slot + 1]);
		s.tracebackDone(aligner);
		drive(slot);
	}
}
