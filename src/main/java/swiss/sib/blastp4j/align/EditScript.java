package swiss.sib.blastp4j.align;

import java.util.Arrays;

/**
 * The operations of a gapped alignment, as runs of one kind ({@code GapEditScript} of BLAST): substitutions, gaps in
 * the query (deletions) and gaps in the subject (insertions).
 */
public final class EditScript {
	/** A gap in the query: consumes subject letters only ({@code eGapAlignDel}). */
	public static final int DEL = 0;
	/** An aligned pair ({@code eGapAlignSub}). */
	public static final int SUB = 3;
	/** A gap in the subject: consumes query letters only ({@code eGapAlignIns}). */
	public static final int INS = 6;

	int[] ops = new int[8];
	int[] nums = new int[8];
	int size;

	public int size() {
		return size;
	}

	public int op(int i) {
		return ops[i];
	}

	public int num(int i) {
		return nums[i];
	}

	/** Appends {@code n} operations, merging with the last run when it has the same kind. */
	void add(int op, int n) {
		if (n == 0) {
			return;
		}
		if (size > 0 && ops[size - 1] == op) {
			nums[size - 1] += n;
			return;
		}
		append(op, n);
	}

	/** Appends a new run without merging. */
	void append(int op, int n) {
		if (size == ops.length) {
			ops = Arrays.copyOf(ops, size * 2);
			nums = Arrays.copyOf(nums, size * 2);
		}
		ops[size] = op;
		nums[size] = n;
		size++;
	}

	void clear() {
		size = 0;
	}

	/** Removes the first run. */
	void removeFirst() {
		System.arraycopy(ops, 1, ops, 0, size - 1);
		System.arraycopy(nums, 1, nums, 0, size - 1);
		size--;
	}

	EditScript copy() {
		EditScript c = new EditScript();
		c.ops = Arrays.copyOf(ops, Math.max(size, 1));
		c.nums = Arrays.copyOf(nums, Math.max(size, 1));
		c.size = size;
		return c;
	}

	/**
	 * {@code Blast_PrelimEditBlockToGapEditScript}: the left extension's operations are already in alignment order,
	 * the right extension's are reversed; the two runs next to the start point are merged when of the same kind.
	 */
	static EditScript join(EditScript left, EditScript right) {
		EditScript e = new EditScript();
		boolean merge = right.size > 0 && left.size > 0 && right.ops[right.size - 1] == left.ops[left.size - 1];
		for (int i = 0; i < left.size; i++) {
			e.append(left.ops[i], left.nums[i]);
		}
		if (right.size == 0) {
			return e;
		}
		int i = right.size - 1;
		if (merge) {
			e.nums[e.size - 1] += right.nums[i];
			i--;
		}
		for (; i >= 0; i--) {
			e.append(right.ops[i], right.nums[i]);
		}
		return e;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < size; i++) {
			sb.append(nums[i]).append(ops[i] == SUB ? 'M' : ops[i] == DEL ? 'D' : 'I');
		}
		return sb.toString();
	}
}
