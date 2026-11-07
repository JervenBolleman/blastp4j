/*
 * qsort with the order of Apple's libc (macOS), for LD_PRELOAD into a Linux BLAST+:
 *
 *     gcc -O2 -shared -fPIC -o apple-qsort.so apple-qsort.c
 *     LD_PRELOAD=$PWD/apple-qsort.so blastp ...
 *
 * Elements that compare equal can end up in a different order than with glibc's qsort, which is a stable merge sort.
 * Apple's qsort is FreeBSD's (Bentley and McIlroy, "Engineering a sort function", 1993) with a recursion depth limit,
 * past which it uses heapsort, and an insertion sort, abandoned after 1 + n / 4 swaps, when a partition made no swaps.
 * This file reimplements that algorithm step by step (from Apple's Libc stdlib/FreeBSD/qsort.c and heapsort.c), so
 * that ties end up in the same order.
 */
#include <stdlib.h>
#include <string.h>

typedef int cmp_t(const void *, const void *);

static void swap(char *a, char *b, size_t es) {
	while (es--) {
		char t = *a;
		*a++ = *b;
		*b++ = t;
	}
}

static char *med3(char *a, char *b, char *c, cmp_t *cmp) {
	return cmp(a, b) < 0 ? (cmp(b, c) < 0 ? b : (cmp(a, c) < 0 ? c : a))
	                     : (cmp(b, c) > 0 ? b : (cmp(a, c) < 0 ? a : c));
}

/* Insertion sort; gives up (returning 0) after more than limit swaps, unless limit is 0. */
static int isort(char *a, size_t n, size_t es, cmp_t *cmp, int limit) {
	int swaps = 0;
	for (char *pm = a + es; pm < a + n * es; pm += es)
		for (char *pl = pm; pl > a && cmp(pl - es, pl) > 0; pl -= es) {
			swap(pl, pl - es, es);
			if (limit && ++swaps > limit)
				return 0;
		}
	return 1;
}

/* Heapsort over elements 1..n: sift down to the bottom, then up (Floyd). */
static void heapsort1(char *v, size_t n, size_t es, cmp_t *cmp) {
	if (n <= 1)
		return;
	char *base = v - es, *k = malloc(es);
	for (size_t l = n / 2; l >= 1; l--)
		for (size_t p = l, c; (c = p * 2) <= n; p = c) {
			if (c < n && cmp(base + c * es, base + (c + 1) * es) < 0)
				c++;
			if (cmp(base + c * es, base + p * es) <= 0)
				break;
			swap(base + p * es, base + c * es, es);
		}
	while (n > 1) {
		memcpy(k, base + n * es, es);
		memcpy(base + n * es, base + es, es);
		n--;
		size_t p = 1, c;
		for (; (c = p * 2) <= n; p = c) {
			if (c < n && cmp(base + c * es, base + (c + 1) * es) < 0)
				c++;
			memcpy(base + p * es, base + c * es, es);
		}
		for (;;) {
			c = p;
			p = c / 2;
			if (c == 1 || cmp(k, base + p * es) < 0) {
				memcpy(base + c * es, k, es);
				break;
			}
			memcpy(base + c * es, base + p * es, es);
		}
	}
	free(k);
}

static void sort(char *a, size_t n, size_t es, cmp_t *cmp, int depth) {
	char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
	size_t d1, d2;
	int r, swapped;
loop:
	if (depth-- <= 0) {
		heapsort1(a, n, es, cmp);
		return;
	}
	if (n <= 7) {
		isort(a, n, es, cmp, 0);
		return;
	}
	pl = a;
	pm = a + (n / 2) * es;
	pn = a + (n - 1) * es;
	if (n > 40) {
		size_t d = (n / 8) * es;
		pl = med3(pl, pl + d, pl + 2 * d, cmp);
		pm = med3(pm - d, pm, pm + d, cmp);
		pn = med3(pn - 2 * d, pn - d, pn, cmp);
	}
	pm = med3(pl, pm, pn, cmp);
	swap(a, pm, es);
	swapped = 0;
	pa = pb = a + es;
	pc = pd = a + (n - 1) * es;
	for (;;) {
		while (pb <= pc && (r = cmp(pb, a)) <= 0) {
			if (r == 0) {
				swapped = 1;
				swap(pa, pb, es);
				pa += es;
			}
			pb += es;
		}
		while (pb <= pc && (r = cmp(pc, a)) >= 0) {
			if (r == 0) {
				swapped = 1;
				swap(pc, pd, es);
				pd -= es;
			}
			pc -= es;
		}
		if (pb > pc)
			break;
		swap(pb, pc, es);
		swapped = 1;
		pb += es;
		pc -= es;
	}
	pn = a + n * es;
	d1 = pa - a < pb - pa ? (size_t)(pa - a) : (size_t)(pb - pa);
	if (d1 > 0)
		for (size_t i = 0; i < d1; i += es)
			swap(a + i, pb - d1 + i, es);
	d1 = (size_t)(pd - pc) < (size_t)(pn - pd - es) ? (size_t)(pd - pc) : (size_t)(pn - pd - es);
	if (d1 > 0)
		for (size_t i = 0; i < d1; i += es)
			swap(pb + i, pn - d1 + i, es);
	if (!swapped && isort(a, n, es, cmp, 1 + n / 4))
		return;
	d1 = pb - pa;
	d2 = pd - pc;
	if (d1 <= d2) {
		if (d1 > es)
			sort(a, d1 / es, es, cmp, depth);
		if (d2 > es) {
			a = pn - d2;
			n = d2 / es;
			goto loop;
		}
	} else {
		if (d2 > es)
			sort(pn - d2, d2 / es, es, cmp, depth);
		if (d1 > es) {
			n = d1 / es;
			goto loop;
		}
	}
}

void qsort(void *a, size_t n, size_t es, cmp_t *cmp) {
	sort(a, n, es, cmp, n ? 2 * (63 - __builtin_clzl(n)) : 0);
}
