package swiss.sib.blastp4j.db;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * A minimal reader of ASN.1 BER, enough for the NCBI binary ASN.1 of BLAST database headers ({@code .phr}): definite
 * and indefinite lengths, constructed context tags, strings and integers.
 */
final class Asn1 {
	/** One element: its tag, and either children (constructed) or content (primitive). */
	record Node(int tagClass, int tag, boolean constructed, List<Node> children, byte[] content) {
		/** The first child with context tag {@code tag}, or null. */
		Node child(int tag) {
			for (Node n : children) {
				if (n.tag == tag && n.tagClass == 2) {
					return n;
				}
			}
			return null;
		}

		/** The single element inside an explicitly tagged element. */
		Node inner() {
			return children.get(0);
		}

		String string() {
			return new String(content, StandardCharsets.ISO_8859_1);
		}

		long integer() {
			long v = content.length > 0 && content[0] < 0 ? -1 : 0;
			for (byte b : content) {
				v = (v << 8) | (b & 0xff);
			}
			return v;
		}

		@Override
		public String toString() {
			if (constructed) {
				return "[" + tagClass + ":" + tag + " " + children + "]";
			}
			return "[" + tagClass + ":" + tag + " '" + string() + "']";
		}
	}

	private Asn1() {
	}

	/** Reads one element at the buffer's position, advancing it. */
	static Node read(ByteBuffer in) {
		int first = in.get() & 0xff;
		int tagClass = first >> 6;
		boolean constructed = (first & 0x20) != 0;
		int tag = first & 0x1f;
		if (tag == 0x1f) { // high tag number
			tag = 0;
			int b;
			do {
				b = in.get() & 0xff;
				tag = (tag << 7) | (b & 0x7f);
			} while ((b & 0x80) != 0);
		}
		int lengthByte = in.get() & 0xff;
		if (lengthByte == 0x80) { // indefinite: children until an end-of-contents marker
			if (!constructed) {
				throw new IllegalArgumentException("indefinite length on a primitive element");
			}
			List<Node> children = new ArrayList<>();
			while (!(in.get(in.position()) == 0 && in.get(in.position() + 1) == 0)) {
				children.add(read(in));
			}
			in.position(in.position() + 2);
			return new Node(tagClass, tag, true, children, null);
		}
		int length = lengthByte;
		if ((lengthByte & 0x80) != 0) {
			length = 0;
			for (int k = lengthByte & 0x7f; k > 0; k--) {
				length = (length << 8) | (in.get() & 0xff);
			}
		}
		if (constructed) {
			int end = in.position() + length;
			List<Node> children = new ArrayList<>();
			while (in.position() < end) {
				children.add(read(in));
			}
			return new Node(tagClass, tag, true, children, null);
		}
		byte[] content = new byte[length];
		in.get(content);
		return new Node(tagClass, tag, false, List.of(), content);
	}
}
