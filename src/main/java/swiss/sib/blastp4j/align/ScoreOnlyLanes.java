package swiss.sib.blastp4j.align;

/**
 * The score-only lanes as their users see them: {@link SemiGappedLanes}, or a copy of it with the gap costs as
 * constants ({@link SemiGappedLanes#of}), a hidden class, which is not a {@link SemiGappedLanes}.
 */
public interface ScoreOnlyLanes {
	/** Extensions submitted and not yet in a lane. */
	int queued();

	/** See {@link SemiGappedLanes#submit(byte[], int, int, byte[], int, int, boolean, int, int)}. */
	void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id);

	/** See {@link SemiGappedLanes#submit(byte[], int, int, byte[], int, int, boolean, int, int, int[])}. */
	void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id,
			int[] matrix);

	/** See {@link SemiGappedLanes#step}. */
	boolean step(SemiGappedLanes.Done done);
}
