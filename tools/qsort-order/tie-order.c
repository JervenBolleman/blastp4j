/*
 * The smallest check that a preloaded qsort orders ties differently, without BLAST+ or a database:
 *
 *     gcc -O2 -shared -fPIC -o apple-qsort.so tools/qsort-order/apple-qsort.c
 *     gcc -O2 -o tie-order tools/qsort-order/tie-order.c
 *     ./tie-order                                 # glibc: 0/0 0/3 0/6 0/9 1/1  1/4 1/7 1/10 2/2  2/5 2/8 2/11
 *     LD_PRELOAD=$PWD/apple-qsort.so ./tie-order  # Apple: 0/0 0/3 0/6 0/9 1/10 1/1 1/7 1/4  2/11 2/5 2/8 2/2
 *
 * Twelve elements with three distinct keys, printed as key/input position. glibc's merge sort leaves equal keys in
 * their input order; Apple's quicksort does not, which is what makes BLAST+'s macOS output differ (see
 * docs/blast-determinism.md). Twelve, because Apple's falls back to an insertion sort at 7 elements or fewer, and
 * that one is stable.
 */
#include <stdio.h>
#include <stdlib.h>

struct item {
	int key, position;
};

static int cmp(const void *a, const void *b) {
	const struct item *x = a, *y = b;
	return (x->key > y->key) - (x->key < y->key);
}

int main(void) {
	struct item v[12];
	for (int i = 0; i < 12; i++) {
		v[i].key = i % 3;
		v[i].position = i;
	}
	qsort(v, 12, sizeof *v, cmp);
	for (int i = 0; i < 12; i++) {
		printf("%d/%d ", v[i].key, v[i].position);
	}
	putchar('\n');
	return 0;
}
