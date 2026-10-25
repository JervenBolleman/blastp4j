package swiss.sib.blastp4j.compo;

import java.util.Arrays;

/**
 * The collection of redone subjects ({@code BlastCompo_Heap} of {@code compo_heap.c}): a plain list until it holds
 * {@code heapThreshold} subjects, then a heap with the worst subject on top. Subjects with an E-value at most
 * {@code ecutoff} are all kept, beyond the threshold; others replace the worst one when better.
 *
 * @param <T> the alignments of a subject
 */
public final class CompoHeap<T> {
	private static final class Rec<T> {
		double bestEvalue;
		int bestScore, subjectIndex;
		T alignments;

		Rec(double bestEvalue, int bestScore, int subjectIndex, T alignments) {
			this.bestEvalue = bestEvalue;
			this.bestScore = bestScore;
			this.subjectIndex = subjectIndex;
			this.alignments = alignments;
		}
	}

	private final int heapThreshold;
	private final double ecutoff;
	private Rec<T>[] array, heapArray; // 1-based
	private int n, capacity;
	private double worstEvalue;

	@SuppressWarnings("unchecked")
	public CompoHeap(int heapThreshold, double ecutoff) {
		this.heapThreshold = heapThreshold;
		this.ecutoff = ecutoff;
		this.capacity = Math.min(100, heapThreshold);
		this.array = new Rec[capacity + 1];
	}

	/** {@code s_CompoHeapRecordCompare}: whether place1 is worse than place2. */
	private static boolean worse(Rec<?> p1, Rec<?> p2) {
		int r = Double.compare(p1.bestEvalue, p2.bestEvalue) > 0 ? 1 : p1.bestEvalue < p2.bestEvalue ? -1 : 0;
		if (r == 0) {
			r = Integer.compare(p2.bestScore, p1.bestScore);
		}
		if (r == 0) {
			r = Integer.compare(p2.subjectIndex, p1.subjectIndex);
		}
		return r > 0;
	}

	private static <T> void swap(Rec<T>[] a, int i, int j) {
		Rec<T> t = a[i];
		a[i] = a[j];
		a[j] = t;
	}

	private static <T> void heapifyDown(Rec<T>[] a, int top, int n) {
		int largest = top, i;
		do {
			i = largest;
			int left = 2 * i, right = 2 * i + 1;
			largest = left <= n && worse(a[left], a[i]) ? left : i;
			if (right <= n && worse(a[right], a[largest])) {
				largest = right;
			}
			if (largest != i) {
				swap(a, i, largest);
			}
		} while (largest != i);
	}

	private static <T> void heapifyUp(Rec<T>[] a, int i) {
		int parent = i / 2;
		while (parent >= 1 && worse(a[i], a[parent])) {
			swap(a, i, parent);
			i = parent;
			parent /= 2;
		}
	}

	private void convertToHeap() {
		if (array != null) {
			heapArray = array;
			array = null;
			for (int i = n / 2; i >= 1; --i) {
				heapifyDown(heapArray, i, n);
			}
		}
	}

	/** {@code BlastCompo_HeapFilledToCutoff}. */
	public boolean filledToCutoff() {
		return n >= heapThreshold && worstEvalue <= ecutoff;
	}

	/** {@code BlastCompo_EarlyTermination} for one query: true when a subject with this E-value can be skipped. */
	public boolean earlyTermination(double evalue) {
		return filledToCutoff() && !(evalue <= 5 * ecutoff);
	}

	/** {@code BlastCompo_HeapWouldInsert}. */
	public boolean wouldInsert(double evalue, int score, int subjectIndex) {
		if (n < heapThreshold || evalue <= ecutoff || evalue < worstEvalue) {
			return true;
		}
		if (heapArray == null) {
			convertToHeap();
		}
		return worse(heapArray[1], new Rec<>(evalue, score, subjectIndex, null));
	}

	private Rec<T>[] insertAtEnd(Rec<T>[] a, Rec<T> r) {
		if (n >= capacity) {
			capacity = Math.max(100 + capacity, (int) (1.5 * capacity));
			a = Arrays.copyOf(a, capacity + 1);
		}
		a[++n] = r;
		return a;
	}

	/** {@code BlastCompo_HeapInsert}. */
	public void insert(T alignments, double evalue, int score, int subjectIndex) {
		if (array != null && n >= heapThreshold) {
			convertToHeap();
		}
		Rec<T> rec = new Rec<>(evalue, score, subjectIndex, alignments);
		if (array != null) {
			array = insertAtEnd(array, rec);
			if (worstEvalue < evalue) {
				worstEvalue = evalue;
			}
		} else {
			if (n < heapThreshold || (evalue <= ecutoff && worstEvalue <= ecutoff)) {
				heapArray = insertAtEnd(heapArray, rec);
				heapifyUp(heapArray, n);
			} else {
				if (worse(heapArray[1], rec)) {
					heapArray[1] = rec;
				}
				heapifyDown(heapArray, 1, n);
			}
			worstEvalue = heapArray[1].bestEvalue;
		}
	}

	/** {@code BlastCompo_HeapPop}: the worst subject's alignments, or null when empty. */
	public T pop() {
		convertToHeap();
		if (n == 0) {
			return null;
		}
		T result = heapArray[1].alignments;
		if (--n > 0) {
			heapArray[1] = heapArray[n + 1];
			heapifyDown(heapArray, 1, n);
		}
		return result;
	}
}
