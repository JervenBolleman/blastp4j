package swiss.sib.blastp4j.util;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;

/**
 * A map of at most {@code capacity} entries that drops the least recently used one when full. Thread safe: the
 * calls are serialized, so keep {@link #computeIfAbsent}'s function short.
 */
public final class LruCache<K, V> {
	private final Map<K, V> map;

	public LruCache(int capacity) {
		if (capacity < 1) {
			throw new IllegalArgumentException("capacity must be at least 1");
		}
		this.map = new LinkedHashMap<>(16, 0.75f, true) {
			@Override
			protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
				return size() > capacity;
			}
		};
	}

	/** The value for {@code key}, made by {@code make} if it is not cached; {@code key} becomes the most recent. */
	public synchronized V computeIfAbsent(K key, Function<? super K, ? extends V> make) {
		return map.computeIfAbsent(key, make);
	}

	public synchronized int size() {
		return map.size();
	}
}
