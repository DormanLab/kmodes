#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <limits.h>

#include "sample.h"
#include "error.h"

void sample(unsigned int N, unsigned int n, unsigned int *idx)
{
	unsigned int t = 0, m = 0;
	double u;

	while (m < n) {
		u = unif_rand();

		if ( (N - t)*u >= n - m )
			++t;
		else
			idx[m++] = t++;
	}
} /* sample */

/**
 * Sample k objects without replacement from n objects with weights.
 *
 * @param n	number of objects
 * @param k	requested sample size
 * @param w	weights
 * @param idx	indices of selected objects
 * @param indic	indicator instead of indices (UINT_MAX if not active)
 * @return	error status
 */
int heap_sample(unsigned int n, unsigned int k, double *w, unsigned int *idx, unsigned int indic)
{
	if (log2(n) > 31)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Sample size exceeded.");

	/* setup heap */
	unsigned int d = 1 << (int)ceil(log2(n));
	double *heap = calloc(2*d, sizeof(*heap));

	if (!heap)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "heap");

	memcpy(&heap[d], w, n*sizeof(*w));
	for (unsigned int i = d - 1; i > 0; --i)
		heap[i] = heap[2*i] + heap[2*i + 1];

	unsigned int offset = d;
	unsigned int j = 0;

	/* sample k objects */
	while (j < k) {
		double r = unif_rand();
		double p = heap[1] * r, left;
		unsigned int i = 1, chosen_idx;

		while (i < offset) {
			i *= 2;
			left = heap[i];
			if (p > left) {
				p -= left;
				i += 1;
			}
		}
		chosen_idx = i - offset;
		if (indic == UINT_MAX) {
			idx[j++] = chosen_idx;
		} else {
			idx[chosen_idx] = indic;
			++j;
		}

		/* reset heap for w[chosen_idx] = 0 */
		i = d + chosen_idx;
		heap[i] = 0;
		while (i > 0) {
			i /= 2;
			heap[i] = heap[2*i] + heap[2*i + 1];
		}
	}

	free(heap);
	return NO_ERROR;
} /* heap_sample */

/**
 * Sample n from N without replacement.  See Vitter, J. S. (1987) ``An Efficient
 * Algorithm for Sequential Random Sampling'' ACM Transactions on Mathematical
 * Software. 13(1): 58--67.
 *
 */
void sample_vitter(unsigned int N, unsigned int n, unsigned int *idx)
{
	if (n >= N)
		return;

	double dn = (double) n;
	double dN = (double) N;
	double ninv = 1.0 / n;
	double vprime = exp(log(unif_rand()) * ninv);
	unsigned int qu1 = 1 + N - n;
	double dqu1 = (double) qu1;
	int nainv = -13.;	/* recommended by Vitter */
	unsigned int thres = -nainv * n;
	unsigned int S=0, limit;
	double nmin1inv, U, X, dnS, y1, y2, top, bottom, V, q;
	unsigned int k = 0;

	while (n > 1 && thres < N) {
		nmin1inv = 1./(-1. + dn);
		do {
			do {
				X = dN * (-vprime + 1.);
				S = (int) X;
				if (S < qu1)
					break;
				vprime = exp(log(unif_rand()) * ninv);
			} while (1);
			U = unif_rand();
			dnS = (double) -S;
			y1 = exp(log(U * dN/dqu1) * nmin1inv);
			vprime = y1 * (-X/dN + 1.) * (dqu1 / (dnS + dqu1));
			if (vprime <= 1.)
				break;
			y2 = 1.;
			top = -1. + dN;
			if (S < n - 1) {
				bottom = -dn + dN;
				limit = N - S;
			} else {
				bottom = dN + dnS - 1.;
				limit = qu1;
			}
			for (unsigned int t = N - 1; t >= limit; --t) {
				y2 = (y2 * top) / bottom;
				top = top - 1.;
				bottom = bottom - 1.;
			}
			if (dN / (dN - X) >= y1 * exp(log(y2)*nmin1inv)) {
				vprime = exp(log(unif_rand()) * nmin1inv);
				break;
			}
			vprime = exp(log(unif_rand()) * ninv);
		} while (1);
		idx[k++] = S;
		N = (N - 1) - S;
		dN = (dN - 1.) + dnS;
		n = n - 1;
		dn = dn - 1.;
		ninv = nmin1inv;
		qu1 = -S + qu1;
		dqu1 = dnS + dqu1;
		thres = thres + nainv;
	}
	if (n > 1) {
		top = N - n;
		while (n >= 2) {
			V = unif_rand();
			q = top / dN;
			while (q > V) {
				S = S + 1;
				top = top - 1.;
				dN = dN - 1.;
				q = (q * top) / dN;
			}
			idx[k++] = S;
			dN = dN - 1.;
			n = n - 1;
		}
	} else {
		S = (int) (N * vprime);
		idx[k++] = S;
	}
} /* sample */

/**
 * Adjust sample_vitter() for ampliclust
 * Sample n from N without replacement.
 *
 * @return		error status
 */
int random_sample(size_t N, size_t n,  size_t *D_idx, size_t *s_idx)
{
	if (n > N)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"invalid random sample");

	size_t t = 0, m = 0;
	double u;

	while (m < n) {
		u = unif_rand();

		if ( (N - t)*u >= n - m )
			++t;
		else
			s_idx[m++] = D_idx[t++];
	}
	return NO_ERROR;
} /* random_sample */
