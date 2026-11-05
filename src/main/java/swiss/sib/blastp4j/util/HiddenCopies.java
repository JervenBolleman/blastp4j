package swiss.sib.blastp4j.util;

import java.io.InputStream;
import java.lang.classfile.ClassFile;
import java.lang.classfile.ClassModel;
import java.lang.classfile.ClassTransform;
import java.lang.classfile.CodeTransform;
import java.lang.classfile.Opcode;
import java.lang.classfile.attribute.InnerClassesAttribute;
import java.lang.classfile.attribute.NestMembersAttribute;
import java.lang.classfile.instruction.FieldInstruction;
import java.lang.constant.ClassDesc;
import java.lang.invoke.MethodHandles;
import java.util.Map;

/**
 * Copies of a class with some of its int (or boolean) fields as constants, defined as hidden classes: the word
 * finders' ({@code seed.WordFinder.of}, {@code seed.BatchWordFinder.of}) and the composition stage's
 * ({@code compo.OptimizeTargetFreq.of}).
 */
public final class HiddenCopies {
	private HiddenCopies() {
	}

	/**
	 * A copy of {@code c}, made with the class-file API and defined as a hidden class in {@code lookup}'s package, in
	 * which each read of a field named in {@code constants} is replaced by its value. The fields are still set, but
	 * no longer read. C2 then folds the values into the compiled code.
	 */
	public static MethodHandles.Lookup define(MethodHandles.Lookup lookup, Class<?> c, Map<String, Integer> constants)
			throws Exception {
		return define(lookup, c, constants, null);
	}

	/**
	 * As {@link #define(MethodHandles.Lookup, Class, Map)}, with {@code more} applied to the code of every method
	 * after the constants (null: none). Create it with {@link CodeTransform#ofStateful} if it keeps state per method.
	 */
	public static MethodHandles.Lookup define(MethodHandles.Lookup lookup, Class<?> c, Map<String, Integer> constants,
			CodeTransform more) throws Exception {
		final byte[] bytes;
		try (InputStream in = c.getResourceAsStream(c.getSimpleName() + ".class")) {
			bytes = in.readAllBytes();
		}
		final ClassFile cf = ClassFile.of();
		final ClassModel model = cf.parse(bytes);
		final ClassDesc self = model.thisClass().asSymbol();
		final CodeTransform replace = (code, element) -> {
			if (element instanceof FieldInstruction f && f.opcode() == Opcode.GETFIELD
					&& f.owner().asSymbol().equals(self) && constants.containsKey(f.name().stringValue())) {
				code.pop(); // the receiver
				code.loadConstant(constants.get(f.name().stringValue()));
			} else {
				code.with(element);
			}
		};
		// A hidden class takes its nest from the lookup, not from these attributes.
		final byte[] copy = cf.transformClass(model,
				ClassTransform.dropping(e -> e instanceof NestMembersAttribute || e instanceof InnerClassesAttribute)
						.andThen(ClassTransform.transformingMethodBodies(more == null ? replace : replace.andThen(more))));
		return lookup.defineHiddenClass(copy, true);
	}
}
