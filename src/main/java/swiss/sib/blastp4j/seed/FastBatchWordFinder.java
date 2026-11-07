package swiss.sib.blastp4j.seed;

import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.util.List;
import java.util.Map;

import swiss.sib.blastp4j.util.HiddenCopies;
import swiss.sib.blastp4j.util.LruCache;

/**
 * The word finder of {@code blastp -task blastp-fast} for a batch of queries: {@link BatchWordFinder} with BLAST's
 * compressed-alphabet scan ({@code s_BlastCompressedAaScanSubject}) and its {@link CompressedTable}. A subject position
 * whose 5 letters are all in the {@link CompressedAlphabet} is looked up; the others (with {@code X}, {@code U},
 * {@code *}, {@code O} or {@code -}) are skipped, as BLAST skips them. The two-hit logic is {@link BatchWordFinder}'s
 * with words of 5 letters, and the hits that are not just stored go to the query's {@link FastWordFinder}.
 * <p>
 * Not thread safe: use one per thread; the {@link CompressedTable} can be shared.
 */
public final class FastBatchWordFinder implements BatchScan {
	private static final int WORD = CompressedTable.WORD;
	/** The weight of a word's newest letter: 15^4. */
	private static final int NEWEST = 50625;
	/** The copies by their constants, as {@link BatchWordFinder}'s. */
	private static final LruCache<Key, MethodHandle> COPIES = new LruCache<>(
			Integer.getInteger("blastp4j.specializations", 256));
	/** Whether this is a GraalVM native image, which cannot define classes at run time. */
	private static final boolean NATIVE_IMAGE = System.getProperty("org.graalvm.nativeimage.imagecode") != null;

	/** A copy's constants. */
	private record Key(int window, int diagMask) {
	}

	private final CompressedTable table;
	private final FastWordFinder[] finders;
	private final int[] diags;
	private final int window, diagMask;
	private int diagOffset;

	public FastBatchWordFinder(byte[][] queries, SeedParameters[] parameters, CompressedTable table) {
		this.table = table;
		window = parameters[0].window();
		diags = new int[table.size()];
		diagMask = table.size() - 1;
		finders = new FastWordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new FastWordFinder(queries[i], parameters[i], diags, table.base()[i]);
		}
		diagOffset = window;
	}

	/**
	 * The word finder for this batch with the window and the diagonal mask as constants, as {@link BatchWordFinder#of}:
	 * a copy of this class shared by the batches with the same values. In a native image, or where the copy fails, a
	 * {@link FastBatchWordFinder}.
	 */
	public static BatchScan of(byte[][] queries, SeedParameters[] parameters, CompressedTable table) {
		if (!NATIVE_IMAGE) {
			try {
				final Key key = new Key(parameters[0].window(), table.size() - 1);
				return (BatchScan) COPIES.computeIfAbsent(key, k -> {
					try {
						final MethodHandles.Lookup lookup = HiddenCopies.define(MethodHandles.lookup(),
								FastBatchWordFinder.class, Map.of("window", k.window(), "diagMask", k.diagMask()));
						return lookup.findConstructor(lookup.lookupClass(), MethodType.methodType(void.class,
								byte[][].class, SeedParameters[].class, CompressedTable.class));
					} catch (Throwable e) {
						throw new IllegalStateException(e);
					}
				}).invoke(queries, parameters, table);
			} catch (IllegalArgumentException e) {
				throw e;
			} catch (Throwable e) {
				// Fall back to the fields.
			}
		}
		return new FastBatchWordFinder(queries, parameters, table);
	}

	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cellStart = table.start(), entry = table.entry(), query = table.query(), base = table.base();
		final long[] present = table.present();
		final int[] d = diags;
		final int diagMask = this.diagMask & (d.length - 1);
		final int offset = diagOffset, windowL = window, lastL = window - WORD - 1;
		final byte[] compress = CompressedAlphabet.COMPRESS;
		final FastWordFinder[] f = finders;
		// The word ending at each position, its oldest letter lowest; good counts the letters since one outside the
		// alphabet.
		int word = 0, good = 0;
		for (int p = 0; p < length; p++) {
			final int c = compress[residues[start + p]];
			if (c >= CompressedAlphabet.SIZE) {
				good = 0;
				word = 0;
				continue;
			}
			word = word / CompressedAlphabet.SIZE + c * NEWEST;
			if (++good < WORD) {
				continue;
			}
			if ((present[word >>> 6] >>> word & 1) == 0) {
				continue; // no hits: start, 3 MB, not read
			}
			final int s = p - (WORD - 1);
			for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
				// As BatchWordFinder.find.
				final int e = entry[k];
				final int diag = (e - s) & diagMask;
				final int state = d[diag];
				final int diff = s - (state - offset), x = diff - WORD;
				if ((~(x | (lastL - x)) | state << 1) >= 0) {
					final int older = (diff - windowL) >> 31;
					d[diag] = (state & older) | ((s + offset) & ~older);
				} else {
					final int i = query[k];
					f[i].hit(residues, start, length, e - base[i], s, out[i]);
				}
			}
		}
		for (FastWordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}
}
