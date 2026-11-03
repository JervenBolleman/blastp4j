package swiss.sib.blastp4j.seed;

import java.io.InputStream;
import java.lang.classfile.ClassFile;
import java.lang.classfile.ClassModel;
import java.lang.classfile.ClassTransform;
import java.lang.classfile.Label;
import java.lang.classfile.MethodModel;
import java.lang.classfile.attribute.InnerClassesAttribute;
import java.lang.classfile.attribute.NestMembersAttribute;
import java.lang.classfile.instruction.SwitchCase;
import java.lang.constant.ClassDesc;
import java.lang.constant.ConstantDescs;
import java.lang.constant.DynamicConstantDesc;
import java.lang.constant.MethodTypeDesc;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.util.ArrayList;
import java.util.List;

/**
 * The word finder for a batch of queries: BLAST scans each subject once for a batch of queries rather than once per
 * query. Each subject position's word is looked up once, in a table that merges the queries' lookup tables, and each
 * hit goes to its query's {@link WordFinder}, with that query's diagonal array, X-drop and cutoff. A query's hits
 * arrive in the order its own word finder takes them (by subject position, then by query offset), so each query gets
 * the alignments {@link WordFinder#find} gives it alone.
 * <p>
 * Not thread safe: use one per thread; the {@link Table} can be shared.
 */
public final class BatchWordFinder implements BatchScan {
	/**
	 * The queries' lookup tables merged: the hits of word {@code w} are entries {@code [start[w], start[w + 1])}, query
	 * {@code query[k]} at offset {@code offset[k]}, by query and then as in that query's table.
	 */
	public record Table(int[] start, int[] query, int[] offset) {
		public static Table of(LookupTable[] tables) {
			int total = 0;
			for (LookupTable t : tables) {
				total += t.offsets().length;
			}
			final int[] start = new int[LookupTable.WORDS + 1], query = new int[total], offset = new int[total];
			int k = 0;
			for (int w = 0; w < LookupTable.WORDS; w++) {
				start[w] = k;
				for (int i = 0; i < tables.length; i++) {
					final int[] ts = tables[i].start(), to = tables[i].offsets();
					for (int j = ts[w]; j < ts[w + 1]; j++) {
						query[k] = i;
						offset[k] = to[j];
						k++;
					}
				}
			}
			start[LookupTable.WORDS] = k;
			return new Table(start, query, offset);
		}
	}

	private final Table table;
	private final Probe[] finders;

	/** The queries, their lookup tables and seeding parameters, by query; {@code table} merges the tables. */
	public BatchWordFinder(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		this.table = table;
		finders = new Probe[queries.length];
		for (int i = 0; i < queries.length; i++) {
			finders[i] = new WordFinder(queries[i], tables[i], parameters[i]);
		}
	}

	private BatchWordFinder(Table table, Probe[] finders) {
		this.table = table;
		this.finders = finders;
	}

	/**
	 * A batch scan whose hits go to copies of {@link WordFinder} from {@link WordFinder#of}, held as constants: a copy
	 * of this class, made with the class-file API and defined as a hidden class, in which {@link #probe} is a switch on
	 * the query to its finder, a constant (class data), so that C2 inlines that finder's {@code probe} with its fields
	 * folded. A copy per call, as the finders keep each scan's diagonals: one per thread. In a native image, or where
	 * the copy fails, a {@link BatchWordFinder}.
	 */
	public static BatchScan of(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		if (!WordFinder.NATIVE_IMAGE) {
			try {
				final Probe[] finders = new Probe[queries.length];
				for (int i = 0; i < queries.length; i++) {
					finders[i] = (Probe) WordFinder.of(queries[i], tables[i], parameters[i]);
				}
				final MethodHandles.Lookup lookup = MethodHandles.lookup()
						.defineHiddenClassWithClassData(specialize(finders.length), List.of(finders), true);
				return (BatchScan) lookup
						.findConstructor(lookup.lookupClass(),
								MethodType.methodType(void.class, Table.class, Probe[].class))
						.invoke(table, finders);
			} catch (Throwable e) {
				// Fall back to the plain finders.
			}
		}
		return new BatchWordFinder(queries, tables, parameters, table);
	}

	/** This class with {@link #probe} a switch on the query to the class data's finder for {@code n} queries. */
	private static byte[] specialize(int n) throws Exception {
		final byte[] bytes;
		try (InputStream in = BatchWordFinder.class.getResourceAsStream("BatchWordFinder.class")) {
			bytes = in.readAllBytes();
		}
		final ClassFile cf = ClassFile.of();
		final ClassModel model = cf.parse(bytes);
		final ClassDesc probe = Probe.class.describeConstable().orElseThrow();
		final ClassDesc list = ConstantDescs.CD_List;
		final MethodTypeDesc probeType = MethodTypeDesc.of(ConstantDescs.CD_void, ConstantDescs.CD_byte.arrayType(),
				ConstantDescs.CD_int, ConstantDescs.CD_int, ConstantDescs.CD_int, ConstantDescs.CD_int, list);
		final ClassTransform switchProbe = (cb, element) -> {
			if (element instanceof MethodModel m && m.methodName().equalsString("probe")) {
				cb.withMethodBody(m.methodName(), m.methodType(), m.flags().flagsMask(), code -> {
					// Locals: finders, query, residues, start, length, q, s, out.
					final Label other = code.newLabel();
					final List<SwitchCase> cases = new ArrayList<>();
					final Label[] labels = new Label[n];
					for (int i = 0; i < n; i++) {
						labels[i] = code.newLabel();
						cases.add(SwitchCase.of(i, labels[i]));
					}
					code.iload(1).tableswitch(0, n - 1, other, cases);
					for (int i = 0; i < n; i++) {
						code.labelBinding(labels[i])
								.ldc(DynamicConstantDesc.ofNamed(ConstantDescs.BSM_CLASS_DATA_AT,
										ConstantDescs.DEFAULT_NAME, probe, i))
								.aload(2).iload(3).iload(4).iload(5).iload(6).aload(7)
								.invokeinterface(probe, "probe", probeType).return_();
					}
					code.labelBinding(other).aload(0).iload(1).aaload().aload(2).iload(3).iload(4).iload(5).iload(6)
							.aload(7).invokeinterface(probe, "probe", probeType).return_();
				});
			} else {
				cb.with(element);
			}
		};
		// A hidden class takes its nest from the lookup, not from these attributes.
		return cf.transformClass(model,
				ClassTransform.dropping(e -> e instanceof NestMembersAttribute || e instanceof InnerClassesAttribute)
						.andThen(switchProbe));
	}

	/** A hit for query {@code i}; in a copy from {@link #of}, a switch to that query's finder. */
	private static void probe(Probe[] finders, int i, byte[] residues, int start, int length, int q, int s,
			List<UngappedHsp> out) {
		finders[i].probe(residues, start, length, q, s, out);
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} for each query, adding
	 * those of query {@code i} to {@code out[i]}, as {@link WordFinder#find} does for one query.
	 */
	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cellStart = table.start(), query = table.query(), offset = table.offset();
		final Probe[] f = finders;
		final int mask = LookupTable.WORDS - 1;
		if (length >= LookupTable.WORD_SIZE) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - LookupTable.WORD_SIZE; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					final int i = query[k];
					probe(f, i, residues, start, length, offset[k], s, out[i]);
				}
			}
		}
		for (Probe finder : f) {
			finder.endSubject(length);
		}
	}
}
