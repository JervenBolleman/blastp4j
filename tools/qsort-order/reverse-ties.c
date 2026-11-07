/*
 * glibc's qsort, then every run of elements that compare equal reversed: the opposite order for every tie. For
 * LD_PRELOAD into a Linux BLAST+, to see how much of its output depends on how ties are ordered:
 *
 *     gcc -O2 -shared -fPIC -o reverse-ties.so reverse-ties.c -ldl
 *     LD_PRELOAD=$PWD/reverse-ties.so blastp ...
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdlib.h>

typedef int cmp_t(const void *, const void *);
typedef void qsort_t(void *, size_t, size_t, cmp_t *);

static void swap(char *a, char *b, size_t es) {
	while (es--) {
		char t = *a;
		*a++ = *b;
		*b++ = t;
	}
}

void qsort(void *v, size_t n, size_t es, cmp_t *cmp) {
	static qsort_t *real;
	if (!real)
		real = (qsort_t *)dlsym(RTLD_NEXT, "qsort");
	real(v, n, es, cmp);
	char *a = v;
	for (size_t i = 0; i < n;) {
		size_t j = i + 1;
		while (j < n && cmp(a + i * es, a + j * es) == 0)
			j++;
		for (size_t l = i, r = j - 1; l < r; l++, r--)
			swap(a + l * es, a + r * es, es);
		i = j;
	}
}
