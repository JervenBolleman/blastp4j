package swiss.sib.blastp4j.seed;

import java.lang.classfile.CodeBuilder;
import java.lang.classfile.CodeElement;
import java.lang.classfile.CodeTransform;
import java.lang.classfile.Instruction;
import java.lang.classfile.Opcode;
import java.lang.classfile.TypeKind;
import java.lang.classfile.instruction.ArrayLoadInstruction;
import java.lang.classfile.instruction.InvokeInstruction;
import java.lang.classfile.instruction.LoadInstruction;
import java.lang.classfile.instruction.StoreInstruction;
import java.lang.constant.ClassDesc;
import java.lang.constant.ConstantDescs;
import java.lang.constant.MethodTypeDesc;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.util.List;
import java.util.Map;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import swiss.sib.blastp4j.util.HiddenCopies;
import swiss.sib.blastp4j.util.LruCache;

/**
 * The word finder for a batch of queries: BLAST scans each subject once for a batch of queries rather than once per
 * query. Each subject position's word is looked up once, in a table that merges the queries' lookup tables, and each
 * hit goes to its query's {@link WordFinder}, with that query's X-drop and cutoff. A query's hits arrive in the order
 * its own word finder takes them (by subject position, then by query offset).
 * <p>
 * The diagonals are BLAST's for the batch: one array for the queries concatenated, as BLAST concatenates them (each
 * query after the previous one and a separator), with the length {@link WordFinder#diagonals} gives for the
 * concatenation. So diagonals of different queries can share an entry, exactly as in BLAST.
 * <p>
 * Not thread safe: use one per thread; the {@link Table} can be shared.
 */
public final class BatchWordFinder implements BatchScan {
	/**
	 * The queries' lookup tables merged: the hits of word {@code w} are entries {@code [start[w], start[w + 1])}, by
	 * query and then as in that query's table. An entry is the offset in the concatenation, {@code base[i]} plus the
	 * query offset, so the diagonal of a hit at subject position {@code s} is {@code (entry - s) & (size - 1)};
	 * {@code query[k]} is the query of entry {@code k}, read only for hits that go on to the two-hit logic.
	 * {@code size} is the length of the diagonal array.
	 * <p>
	 * {@code narrow} holds the same entries as chars, half the memory the scan reads, when the concatenation is
	 * shorter than 65,536 letters (BLAST's batches are about 10,000), and is null otherwise: the scan's copy for such a
	 * table reads it ({@link BatchWordFinder#of}).
	 */
	public record Table(int[] start, int[] entry, char[] narrow, int[] query, int[] base, int size) {
		public static Table of(byte[][] queries, LookupTable[] tables, int window) {
			final int n = queries.length;
			// BLAST's context offsets (s_QueryInfo_SetContext): one separator after each query that is not empty.
			final int[] base = new int[n];
			long length = 0;
			for (int i = 0; i < n; i++) {
				if (i > 0 && queries[i - 1].length > 0) {
					length++;
				}
				base[i] = (int) length;
				length += queries[i].length;
			}
			if (length + window > 1 << 30) {
				throw new IllegalArgumentException("a batch of " + n + " queries is too long: " + length);
			}
			int total = 0;
			for (LookupTable t : tables) {
				total += t.offsets().length;
			}
			// The entries are padded with a vector of zeros, so that the scan's lanes can load past a word's last one.
			final int[] start = new int[LookupTable.WORDS + 1], entry = new int[total + LANES], query = new int[total];
			int k = 0;
			for (int w = 0; w < LookupTable.WORDS; w++) {
				start[w] = k;
				for (int i = 0; i < n; i++) {
					final int[] ts = tables[i].start(), to = tables[i].offsets();
					for (int j = ts[w]; j < ts[w + 1]; j++) {
						query[k] = i;
						entry[k++] = base[i] + to[j];
					}
				}
			}
			start[LookupTable.WORDS] = k;
			char[] narrow = null;
			if (length <= Character.MAX_VALUE) {
				narrow = new char[total];
				for (int j = 0; j < total; j++) {
					narrow[j] = (char) entry[j];
				}
			}
			return new Table(start, entry, narrow, query, base, WordFinder.diagonals((int) length, window));
		}
	}

	/**
	 * The copies by their constants. Nearly every full batch has the same ones: BLAST's batches are about 10,000
	 * letters, so their diagonal mask is 16383.
	 */
	private static final LruCache<Key, MethodHandle> COPIES = new LruCache<>(
			Integer.getInteger("blastp4j.specializations", 256));
	/** Whether this is a GraalVM native image, which cannot define classes at run time. */
	private static final boolean NATIVE_IMAGE = System.getProperty("org.graalvm.nativeimage.imagecode") != null;
	private static final VectorSpecies<Integer> INTS = IntVector.SPECIES_PREFERRED;
	private static final int LANES = INTS.length();
	/**
	 * Whether the scan takes a subject position's hits a vector at a time (needs 8 int lanes with a gather, as
	 * {@code SemiGappedLanes}). The lanes then read the int entries: the copy does not read the narrow ones.
	 */
	public static final boolean SEED_LANES = LANES >= 8 && !NATIVE_IMAGE
			&& Boolean.parseBoolean(System.getProperty("blastp4j.seedLanes", "true"));

	/** A copy's constants, and whether it reads the table's narrow entries. */
	private record Key(int window, int diagMask, boolean narrow) {
	}

	private final Table table;
	private final WordFinder[] finders;
	/** The diagonals of the concatenated queries. */
	private final int[] diags;
	private final int window, diagMask;
	/** As each finder's {@code diagOffset}, which all move together. */
	private int diagOffset;
	/** The lanes' diagonals and their new states, stored to {@link #diags} one lane at a time. */
	private final int[] laneDiag = new int[LANES], laneState = new int[LANES];

	/** The queries, their lookup tables and seeding parameters, by query; {@code table} merges the tables. */
	public BatchWordFinder(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		this.table = table;
		window = parameters[0].window();
		diags = new int[table.size()];
		diagMask = table.size() - 1;
		finders = new WordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new WordFinder(queries[i], tables[i], parameters[i], diags, table.base()[i]);
		}
		diagOffset = window;
	}

	/**
	 * The word finder for this batch, with the window and the diagonal mask as constants: a copy of this class (see
	 * {@link WordFinder#of}), shared by the batches with the same values. Only the scan is copied: the queries'
	 * finders are {@link WordFinder}s, so that all hits go to the one class. Where the table has narrow entries, the
	 * copy's scan reads those instead ({@link #narrowEntries}). In a native image, or where the copy fails, a
	 * {@link BatchWordFinder}, which reads the int entries.
	 */
	public static BatchScan of(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		if (!NATIVE_IMAGE) {
			try {
				final Key key = new Key(parameters[0].window(), table.size() - 1, table.narrow() != null && !SEED_LANES);
				return (BatchScan) COPIES.computeIfAbsent(key, k -> {
					try {
						final MethodHandles.Lookup lookup = HiddenCopies.define(MethodHandles.lookup(),
								BatchWordFinder.class, Map.of("window", k.window(), "diagMask", k.diagMask()),
								k.narrow() ? narrowEntries() : null);
						return lookup.findConstructor(lookup.lookupClass(), MethodType.methodType(void.class,
								byte[][].class, LookupTable[].class, SeedParameters[].class, Table.class));
					} catch (Throwable e) {
						throw new IllegalStateException(e);
					}
				}).invoke(queries, tables, parameters, table);
			} catch (IllegalArgumentException e) {
				throw e;
			} catch (Throwable e) {
				// Fall back to the fields.
			}
		}
		return new BatchWordFinder(queries, tables, parameters, table);
	}

	/**
	 * The code transform of {@link #of}'s copy for a table with narrow entries: in {@link #find}, {@code table.entry()}
	 * becomes {@code table.narrow()}, and each load from it ({@code aload entry; iload k; iaload}, as javac compiles
	 * {@code entry[k]}) a char load. A load of the array that is not followed so leaves an int load of a char array,
	 * which the verifier rejects: then {@link #of} falls back to this class.
	 */
	private static CodeTransform narrowEntries() {
		final ClassDesc table = ClassDesc.of(Table.class.getName());
		final MethodTypeDesc entryType = MethodTypeDesc.of(ConstantDescs.CD_int.arrayType());
		final MethodTypeDesc narrowType = MethodTypeDesc.of(ConstantDescs.CD_char.arrayType());
		return CodeTransform.ofStateful(() -> new CodeTransform() {
			/** The local of the entry array (-1: none yet), whether its store is next, and the load pattern's step. */
			private int slot = -1, step;
			private boolean store;

			@Override
			public void accept(CodeBuilder code, CodeElement e) {
				if (e instanceof InvokeInstruction i && i.opcode() == Opcode.INVOKEVIRTUAL
						&& i.owner().asSymbol().equals(table) && i.name().equalsString("entry")
						&& i.typeSymbol().equals(entryType)) {
					code.invokevirtual(table, "narrow", narrowType);
					store = true;
					return;
				}
				if (store && e instanceof StoreInstruction st) {
					slot = st.slot();
					store = false;
				} else if (e instanceof LoadInstruction l && l.typeKind() == TypeKind.REFERENCE && l.slot() == slot) {
					step = 1;
				} else if (step == 1 && e instanceof LoadInstruction l && l.typeKind() == TypeKind.INT) {
					step = 2;
				} else if (step == 2 && e instanceof ArrayLoadInstruction a && a.typeKind() == TypeKind.INT) {
					code.arrayLoad(TypeKind.CHAR);
					step = 0;
					return;
				} else if (e instanceof Instruction) {
					step = 0;
				}
				code.with(e);
			}
		});
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} for each query, adding
	 * those of query {@code i} to {@code out[i]}, as {@link WordFinder#find} does for one query.
	 */
	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		if (SEED_LANES) {
			findLanes(residues, start, length, out);
			return;
		}
		final int[] cellStart = table.start(), entry = table.entry(), query = table.query(), base = table.base();
		final int[] d = diags;
		final int diagMask = this.diagMask;
		final int offset = diagOffset, windowL = window, lastL = window - LookupTable.WORD_SIZE - 1;
		final WordFinder[] f = finders;
		final int mask = LookupTable.WORDS - 1;
		if (length >= LookupTable.WORD_SIZE) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - LookupTable.WORD_SIZE; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					// The usual hits, decided from the entry: the diagonal's last hit is a window or more back, so
					// this one replaces it, or it overlaps this one, which is then ignored. Only hits to extend and
					// hits on an extended diagonal go to the query's finder.
					final int e = entry[k];
					// Also masked with the array's length (a power of two), so that C2 knows the index is in bounds.
					final int diag = (e - s) & diagMask & (d.length - 1);
					final int state = d[diag];
					// Without branches but one: x = diff - WORD_SIZE is in [0, last] (a hit to extend) when neither
					// x nor last - x is negative, and the extended flag (bit 30; bit 31 is never set) goes to the
					// sign bit.
					final int diff = s - (state - offset), x = diff - LookupTable.WORD_SIZE;
					if ((~(x | (lastL - x)) | state << 1) >= 0) {
						final int older = (diff - windowL) >> 31; // -1 if the hit overlaps the last one: keep that
						d[diag] = (state & older) | ((s + offset) & ~older);
					} else {
						final int i = query[k];
						f[i].hit(residues, start, length, e - base[i], s, out[i]);
					}
				}
			}
		}
		for (WordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}

	/**
	 * {@link #find} with a subject position's hits a vector at a time. The hits of one position are on different
	 * diagonals (the diagonal array is longer than the concatenation), so a vector of them can be decided at once: the
	 * same decisions as {@link #find}, and the hits that go to a query's finder go in the same order.
	 */
	private void findLanes(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cellStart = table.start(), entry = table.entry(), query = table.query(), base = table.base();
		final int[] d = diags, ld = laneDiag, ls = laneState;
		final int diagMask = this.diagMask & (d.length - 1);
		final int offset = diagOffset, windowL = window, lastL = window - LookupTable.WORD_SIZE - 1;
		final WordFinder[] f = finders;
		final int mask = LookupTable.WORDS - 1;
		final IntVector iota = IntVector.zero(INTS).addIndex(1);
		if (length >= LookupTable.WORD_SIZE) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - LookupTable.WORD_SIZE; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				final int end = cellStart[word + 1];
				for (int k = cellStart[word]; k < end; k += LANES) {
					final VectorMask<Integer> valid = iota.compare(VectorOperators.LT, end - k);
					final IntVector diag = IntVector.fromArray(INTS, entry, k).sub(s).and(diagMask);
					diag.intoArray(ld, 0);
					final IntVector state = IntVector.fromArray(INTS, d, 0, ld, 0);
					// As in find: x in [0, last] and not extended is a hit to extend; the others are usual hits.
					final IntVector diff = state.neg().add(s + offset), x = diff.sub(LookupTable.WORD_SIZE);
					final VectorMask<Integer> usual = x.or(x.neg().add(lastL)).not().or(state.lanewise(
							VectorOperators.LSHL, 1)).compare(VectorOperators.GE, 0).and(valid);
					final IntVector older = diff.sub(windowL).lanewise(VectorOperators.ASHR, 31);
					final IntVector replaced = state.and(older).or(older.not().and(s + offset));
					// The other lanes write back what they read; stored from the last lane down, so that a padding
					// lane on the diagonal of a hit never overwrites that hit's new state.
					state.blend(replaced, usual).intoArray(ls, 0);
					for (int l = LANES - 1; l >= 0; l--) {
						d[ld[l] & diagMask] = ls[l];
					}
					long rest = valid.toLong() & ~usual.toLong();
					for (; rest != 0; rest &= rest - 1) {
						final int j = k + Long.numberOfTrailingZeros(rest), i = query[j];
						f[i].hit(residues, start, length, entry[j] - base[i], s, out[i]);
					}
				}
			}
		}
		for (WordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}
}
