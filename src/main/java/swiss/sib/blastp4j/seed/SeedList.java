package swiss.sib.blastp4j.seed;

import java.util.AbstractList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Objects;
import java.util.RandomAccess;

/**
 * The ungapped alignments the word finders find, a list over an array that grows by {@link #GROW} elements when full.
 * The finders append to {@link #items} directly ({@code WordFinder.hit}), so that saving an alignment is no call but
 * on the rare path where the array is full: C2 does not inline {@code ArrayList.add} where it is rarely reached, and
 * a call keeps the values of the loop around it on the stack (C2 has no callee-saved registers).
 */
public final class SeedList extends AbstractList<UngappedHsp> implements RandomAccess {
	/** Elements added when the array is full. */
	static final int GROW = 64;

	UngappedHsp[] items = new UngappedHsp[GROW];
	int size;

	/** Makes room for {@link #GROW} more elements. */
	void grow() {
		items = Arrays.copyOf(items, items.length + GROW);
	}

	@Override
	public UngappedHsp get(int index) {
		Objects.checkIndex(index, size);
		return items[index];
	}

	@Override
	public UngappedHsp set(int index, UngappedHsp element) {
		Objects.checkIndex(index, size);
		final UngappedHsp old = items[index];
		items[index] = element;
		return old;
	}

	@Override
	public int size() {
		return size;
	}

	@Override
	public boolean add(UngappedHsp h) {
		if (size == items.length) {
			grow();
		}
		items[size++] = h;
		return true;
	}

	@Override
	public void clear() {
		Arrays.fill(items, 0, size, null);
		size = 0;
	}

	/** Stable, as {@code ArrayList.sort}. */
	@Override
	public void sort(Comparator<? super UngappedHsp> c) {
		Arrays.sort(items, 0, size, c);
	}
}
