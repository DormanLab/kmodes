//#define DEBUGGING_CODE

#include <Rmath.h>
#include <limits.h>
#ifdef DEBUGGING_CODE
#include <stdio.h>
#endif

#include "entropy.h"
#include "run_kmodes.h"
#include "error.h"
#include "order.h"
#include "hash.h"
#include "sample.h"

static inline double hd(data_t *x, data_t *y, unsigned int p);
static int perturb_by_hd(data *dat, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
static int initialize_proportional_to_abundance(data *dat, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
static double entropy(data_t *x, int unique_c, int len);
static int entropy_per_coordinate(data *dat);
static int mask_idx(double *entr, unsigned int p, double percentage, size_t **indices);
unsigned int binary_search(data_t **x, data *dat, size_t *indices, int l, int r, unsigned int K, int abunk);
int randomize_masked_hash(data *dat);
int select_by_hd (data *dat, data_t **seeds, data_t **nseeds, unsigned int *nsd_idx, unsigned int *sd_idx, unsigned int K);
void replace_by_random(data *dat, unsigned int K, unsigned int k1, data_t **seeds,
					   unsigned int *sd_idx);

static inline double hd(data_t *x, data_t *y, unsigned int p)
{
	double d = 0;
	for (unsigned int j = 0; j < p; ++j)
		d += (x[j] != y[j]);
	return d;
} /* hd */

/**
 * Compute entropy of categorical data.
 *
 * @param x		categorical data (1 x n)
 * @param unique_c	number of categories
 * @param n		number of observations, n
 * @return		computed entropy
 */
static double entropy(data_t *x, int unique_c, int n)
{
	double entro = 0;
	int hist[unique_c];

	/* assume categories are 0, 1, 2, ... without skips */
	for (int i = 0; i < unique_c; ++i)
		hist[i] = 0;
	for (int i = 0; i < n; ++i)
		hist[x[i]]++;
	
	for (int i = 0; i < unique_c; ++i) {
		if (hist[i] == 0)
			continue;
		entro -= hist[i] * (log2(hist[i]) - log2(n)) / n;
	}
	
	return entro;
} /* entropy */

/**
 * Compute entropy of each coordinate.
 *
 * @param dat		data object
 * @return		error status
 */
static int entropy_per_coordinate(data *dat)
{
	unsigned int n = dat->n_observations, p = dat->n_coordinates;
	unsigned int i, j;
	data_t max_ncat = dat->max_n_categories;
	data_t **x = dat->dmat;

	dat->entropy = malloc(p * sizeof(*dat->entropy));
	if (!dat->entropy)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data::entropy");
	dat->total_entropy = 0;

#ifdef DEBUGGING_CODE
	fprintf(stderr, "Entropy: ");
#endif
	for (j = 0; j < p; ++j) {
		data_t tmp[n];

		for (i = 0; i < n; ++i)
			tmp[i] = x[i][j];

		dat->entropy[j] = entropy(tmp, max_ncat, n);
		dat->total_entropy += dat->entropy[j];
#ifdef DEBUGGING_CODE
		fprintf(stderr, "%lf ", dat->entropy[j]);
#endif
	}
#ifdef DEBUGGING_CODE
	fprintf(stderr, "\n");
#endif

	return NO_ERROR;
} /* entropy_per_coordinate */

/**
 * Order of coordinates in increasing entropy order.
 *
 * @param s_entropy	entropy across coordinates
 * @param p		number of coordinates
 * @param percentage	not used
 * @param indices	order of coordinates sorted by entropy
 * @return		error status
 */
static int mask_idx(double *s_entropy, unsigned int p, double percentage, size_t **indices)
{
//	size_t *pos = NULL;
	UNUSED(percentage);
//	*indices = malloc(sizeof(size_t) * p);
	
	*indices = qorder_double(s_entropy, p); // get the order of entr

	/*
	for (unsigned int j = 0, i = p - 1; j < i; ++j, --i) {
		size_t tmp = (*indices)[j];
		(*indices)[j] = (*indices)[i];
		(*indices)[i] = tmp;
	}
	*/
#ifdef DEBUGGING_CODE
	fprintf(stderr, "Ranked sites: ");
	for (unsigned int j = 0; j < p; ++j)
		fprintf(stderr, "%zu ", (*indices)[j]);
	fprintf(stderr, "\n");
#endif
//	free(pos);

	return NO_ERROR;
} /* mask_idx */

/**
 * Binary search for masked sites.
 *
 * @param x		data (n x p)
 * @param dat		data object
 * @param indices	??
 * @param l		left bounding index
 * @param r		right bounding index
 * @param K		number of clusters
 * @param abunk		desired number of masked abundances >1
 * @return		entropy order statistic that yields
 */
unsigned int binary_search(data_t **x, data *dat, size_t *indices,
				 int l, int r, unsigned int K, int abunk)
{
	unsigned int i, j;
	unsigned int n = dat->n_observations, p = dat->n_coordinates;
	unsigned int mid = 0;

	if (r > l) {
		mid = l + (r - l) / 2;
		// this is not efficient
		for (j = l; j < mid; ++j)
			for (i = 0; i < n; ++i)
				x[i][indices[j]] = 0;
		
		/* build hash table */
		dat->hash_length = 0;
		dat->seq_count = NULL;	/* [KSD,BUG] was memory leak */
		for (i = 0; i < n; ++i)
			dat->hash_length
				+= add_sequence(&dat->seq_count, x[i], p, i);
		int pool = 0;
		for (hash *s = dat->seq_count; s != NULL; s = s->hh.next)
			if (s->count > 1)
				pool++;
		if (pool == abunk) {
			return mid;
		} else if (pool < abunk) {// mask more
			delete_all(dat->seq_count);
			for (i = 0; i < n; ++i)
				x[i][indices[mid]] = 0;
			return binary_search(x, dat, indices, mid + 1, r, K,
									abunk);
		} else if (l < mid) {	// mask less, revert copied first
			for (j = l; j < mid; ++j)
				for (i = 0; i < n; ++i)
					x[i][indices[j]]
						= dat->dmat[i][indices[j]];
			delete_all(dat->seq_count);
			return binary_search(x, dat, indices, l, mid, K, abunk);
		} else {
			return mid;
		}
	} else {	/* r == l */
		dat->hash_length = 0;
		dat->seq_count = NULL;
		for (i = 0; i < n; ++i)
			dat->hash_length
				+= add_sequence(&dat->seq_count, x[i], p, i);
		return l;
	}
} /* binary_search */

/**
 * Create as masked hash.
 *
 * @param dat	data object
 * @return	error status
 */
int mask_nhash(data *dat, options *opt)
{
	int err = NO_ERROR;
	unsigned int K = opt->K;
	int abunk = opt->abunk;
	unsigned int p = dat->n_coordinates, n = dat->n_observations;
	unsigned int i, m;
	
	if ((err = entropy_per_coordinate(dat)))
		return err;

	err = mask_idx(dat->entropy, p, 1, &dat->entropy_order);

	/* [KSD,TODO] Crazy expensive to copy data! */
	dat->masked_dmat = malloc(n * sizeof(*dat->masked_dmat));

	if (!dat->masked_dmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"data::masked_dmat");

	dat->masked_dmat[0] = malloc(n * p * sizeof(**dat->masked_dmat));

	if (!dat->masked_dmat[0])
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"data::masked_dmat[0]");
	
	memcpy(dat->masked_dmat[0], dat->data, n * p * sizeof(*dat->data));
	for (i = 1; i < n; ++i)
		dat->masked_dmat[i] = dat->masked_dmat[0] + i * p;

	// use binary search to find out the number of masked position
	dat->n_masked = binary_search(dat->masked_dmat, dat,
					dat->entropy_order, 0, p - 1, K, abunk);
#ifdef DEBUGGING_CODE
	fprintf(stderr, "Masked sites (%u): ", dat->n_masked);
	for (j = 0; j < dat->n_masked; ++j)
		fprintf(stderr, "%zu ", dat->entropy_order[j]);
	fprintf(stderr, "\n");
#endif

	/* store index of reads for all unique sequences */
	for (i = 0; i < n; ++i)
		if ((err = add_seq_idx(dat->seq_count,
						dat->masked_dmat[i], p, i)))
			return err;

	sort_by_count(&dat->seq_count);

	m = 0;
	for (hash *s = dat->seq_count; s != NULL; s = s->hh.next, ++m) {
#ifdef DEBUGGING_CODE
		fprintf(stderr, "Abundance[%zu(%u)=", s->idx, m);
		for (unsigned int j = 0; j < p; ++j)
			if (dat->entropy_order[j] < dat->n_masked)
				fprintf(stderr, "[%d]", dat->masked_dmat[s->idx][j]);
			else
				fprintf(stderr, "%d", dat->masked_dmat[s->idx][j]);
#endif
		if (opt->true_modes) {
			unsigned int min_hd = UINT_MAX;
#ifdef DEBUGGING_CODE
			unsigned int chosen_k = 0;
#endif

			for (unsigned int k = 0; k < opt->true_K; ++k) {
				unsigned int hd = 0;

				for (unsigned int j = 0; j < p; ++j)
					if (dat->entropy_order[j] >= dat->n_masked
						&& dat->dmat[s->idx][j] != opt->true_modes[k][j])
						++hd;
				/*
				if (!hd) {
					fprintf(stderr, " matches true mode %u", k);
					break;
				} else*/ if (hd < min_hd) {
					min_hd = hd;
#ifdef DEBUGGING_CODE
					chosen_k = k;
#endif
				}
			}
#ifdef DEBUGGING_CODE
			fprintf(stderr, " HD=%u%s from true mode %u", min_hd, !min_hd?"***":"", chosen_k);
#endif
		}
#ifdef DEBUGGING_CODE
		fprintf(stderr, "]: %u\n", s->count);
#endif
	}
	
	dat->mask = calloc(dat->n_coordinates, sizeof(*dat->mask));

	if (!dat->mask)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "mask");

	return err;
} /* mask_nhash */

/**
 * Select sites to leave unmasked in proportion to their entropy.
 *
 * @param dat	data object
 * @return	error status
 */
int randomize_masked_hash(data *dat)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
#endif
	int err;
	unsigned int n = dat->n_observations, p = dat->n_coordinates, i, j;
	unsigned int n_unmasked = p - dat->n_masked;

	/* reset */
#ifdef DEBUGGING_CODE
	debug_msg(fxn_debug >= DEBUG_I, DEBUG_I, "Reseting mask: will mask %u of %u.\n", dat->n_masked, p);
#endif
	delete_all(dat->seq_count);
	if (n_unmasked == p)
		memset(dat->mask, 0, p * sizeof(*dat->mask));
	else
		memset(dat->mask, -1, p * sizeof(*dat->mask));
	memcpy(dat->masked_dmat[0], dat->data, n * p * sizeof(*dat->data));

#ifdef DEBUGGING_CODE
	debug_msg(fxn_debug >= DEBUG_I, DEBUG_I, "Mask:");
	for (j = 0; j < p; ++j)
		debug_msg_cont(fxn_debug >= DEBUG_I, DEBUG_I, " %zu=%f (%d)", dat->entropy_order[j], dat->entropy[dat->entropy_order[j]], dat->mask[dat->entropy_order[j]]);
	debug_msg_cont(fxn_debug >= DEBUG_I, DEBUG_I, "\n");
#endif

	if (n_unmasked < p) {
		/* choose masked sites proportionally to entropy */
		heap_sample(p, n_unmasked, dat->entropy, dat->mask, 0);

#ifdef DEBUGGING_CODE
		fprintf(stderr, "Masked sites (%u):", dat->n_masked);
#endif

		/* mask the masked sites */
		for (j = 0; j < p; ++j) {
			if (!dat->mask[j])
				continue;
#ifdef DEBUGGING_CODE
			fprintf(stderr, " %u", j);
#endif
			/* mask chosen position */
			for (i = 0; i < n; ++i)
				dat->masked_dmat[i][j] = 0;
		}
#ifdef DEBUGGING_CODE
		fprintf(stderr, "\n");
#endif
	}

	/* create the masked hashed */
	dat->hash_length = 0;
	dat->seq_count = NULL;
	for (i = 0; i < n; ++i)
		dat->hash_length += add_sequence(&dat->seq_count,
						dat->masked_dmat[i], p, i);
#ifdef DEBUGGING_CODE
		debug_msg(fxn_debug >= DEBUG_I, DEBUG_I, "Created hash of length %u.\n", dat->hash_length);
#endif
	for (i = 0; i < n; ++i)
		if ((err = add_seq_idx(dat->seq_count,
						dat->masked_dmat[i], p, i)))
			return err;

	return NO_ERROR;
} /* randomize_masked_hash */


/**
 * Init with abundance by entropy mask high entropy site
 *
 * @param dat		data object
 * @param K		number of clusters
 * @param k1		number of fixed seeds
 * @param seeds		the seeds (to set)
 * @param sd_idx	seed indices
 * @param inner_iter	inner iteration number
 * @return		error status
 */
int
initialize_by_abundance(
				data *dat,
				unsigned int K,
				unsigned int k1,
				data_t **seeds,
				unsigned int *sd_idx)
{
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
#endif
	int err = NO_ERROR;
	unsigned int p = dat->n_coordinates;
	unsigned int j, k, l, nties = 0;
	data_t **x = dat->dmat;
	hash *sc = NULL, *snext = NULL;

	if ((err = randomize_masked_hash(dat)))
		return err;
	
	sort_by_count(&dat->seq_count);
	sc = dat->seq_count;

	for (j = k1; j < K;) {

		/* consider next most abundant seed */
		//unsigned int s = sc->idx;
	   
		/* count ties */
		snext = sc->hh.next;
		if (!snext)
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "There are "
				"fewer than %u distinct sequences.\n", K);
		nties = 1;
		while (snext->count == sc->count) {
			++nties;
			snext = snext->hh.next;
			if (!snext)
				break;
		}

		/* check which of these are already selected */
		unsigned char used[nties];
		unsigned int nused = 0;

		snext = sc;
		for (l = 0; l < nties; ++l) {
			used[l] = 0;
			for (k = 0; k < k1; ++k)
				if (!hd(x[k], x[snext->idx], p)) {
					++nused;
					used[l] = 1;
					break;
				}
			snext = snext->hh.next;
		}

#ifdef DEBUGGING_CODE
		debug_msg(DEBUG_I <= fxn_debug, DEBUG_I, "Abundance %u, %u ties (%u used)\n", sc->count, nties, nused);
#endif

		/* all are used */
		if (nties == nused) {
			for (l = 0; sc != NULL && l < nties; ++l)
				sc = sc->hh.next;
			continue;
		}

		/* can use them all */
		if (K - j >= nties - nused) {
			for (l = 0; l < nties; ++l) {
				if (!used[l]) {
#ifdef DEBUGGING_CODE
					if (DEBUG_II <= fxn_debug) {
						for (unsigned int m = 0; m < p; ++m)
							fprintf(stderr, "%d", x[sc->idx][m]);
						fprintf(stderr, "\n");
					}
#endif
					/* avoid identical seeds */
					for (k = 0; k < j; ++k)
						if (!hd(seeds[k], x[sc->idx], p * sizeof(**seeds)))
							break;
					if (k < j)
						mmessage(ERROR_MSG, INTERNAL_ERROR, "Identical seeds!\n");
					//randomly sample one from the selected pool
					int r = (int) runif(0, sc->count);
					size_t id = sc->idx_array[r];

					memcpy(seeds[j], x[id], p * sizeof(**seeds));
					if (sd_idx)
						sd_idx[j] = sc->idx;
					++j;
				}
				sc = sc->hh.next;
			}
			continue;
		}

		/* need to randomly sample those to keep */
		while (j < K && nties - nused) {
			/* randomly select an unused one */
			k = (unsigned int) runif(0, nties);

			if (used[k])
				continue;

			/* forward to it */
			snext = sc;
			for (l = 0; l < k; ++l)
				snext = snext->hh.next;
			/* avoid identical members */
			for (l = 0; l < j; ++l)
				if (!hd(seeds[l], x[snext->idx], p * sizeof(**seeds)))
					break;
			if (l < j)
				mmessage(ERROR_MSG, INTERNAL_ERROR, "Identical seeds!\n");
#ifdef DEBUGGING_CODE
			if (DEBUG_II <= fxn_debug) {
				for (unsigned int m = 0; m < p; ++m)
					fprintf(stderr, "%d", x[snext->idx][m]);
				fprintf(stderr, "\n");
			}
#endif
			int r = runif(0, snext->count);
			size_t id = snext->idx_array[r];

			memcpy(seeds[j], x[id], p * sizeof(**seeds));

			if (sd_idx)
				sd_idx[j] = snext->idx;
			++j;

			/* record it as used */
			used[k] = 1;
			++nused;
		}

		for (l = 0; sc != NULL && l < nties; ++l)
			sc = sc->hh.next;
	}
	
	return err;
} /* initialize_by_abundance */

/**
 * Perturb for by masked abundance.
 *
 * @param dat	data object
 * @param opt	options object
 * @return	error status
 */
int perturb(data *dat, options *opt)
{
	return perturb_by_hd(dat, opt->K,
		opt->n_seed_set, dat->seeds, dat->seed_idx);
} /* perturb */

/**
 * Initialize proportional to abundance. Same as random sampling
 * except with k1 > 0.
 */
static int
initialize_proportional_to_abundance(
				data *dat,
				unsigned int K,
				unsigned int k1,
				data_t **seeds,
				unsigned int *sd_idx)
{
	data_t **x = dat->dmat;
	unsigned int p = dat->n_coordinates;
	double dsum = 0, r;
	unsigned char repeat = 0;
	unsigned int s = 0;
	hash *useq;

	for (useq = dat->seq_count; useq != NULL; useq = useq->hh.next)
		dsum += useq->count;

	for (unsigned int k = k1; k < K; ++k) {
		if (repeat)
			--k;
		repeat = 0;
		r = runif(0, dsum);

		for (useq = dat->seq_count, dsum = useq->count; useq != NULL && r > dsum; useq = useq->hh.next, dsum += useq ? useq->count : 0);

		s = useq->idx;
		// get the center of obs with the chosen masked obs
		int r = runif(0, useq->count);
		size_t id = useq->idx_array[r];

		memcpy(seeds[k], x[id], p * sizeof(**seeds));
		sd_idx[k] = s;

		for (unsigned int l = 0; l < k1; ++l)
			if (!hd(seeds[l], x[s], p)) {
				repeat = 1;
				break;
			}
		
//		s = dat->cluster_id[s];
//		memcpy(seeds[k], dat->seeds[s], p * sizeof(**seeds));
//		sd_idx[k] = s;
//
//		for (unsigned int l = 0; l < k1; ++l)
//			if (!hd(seeds[l], dat->seeds[s], p)) {
//				repeat = 1;
//				break;
//			}
	}

	return NO_ERROR;
} /* initialize_proportional_to_abundance */

/**
 *
 * @param dat		data object (k-modes)
 * @param K		number of clusters
 * @param k1		unused
 * @param seeds		the current seeds
 * @param sd_idx	the current seed indices
 * @return		error status
 */
static int
perturb_by_hd(
				data *dat,
				unsigned int K,
				unsigned int k1,
				data_t **seeds,
				unsigned int *sd_idx)
{																		   /**/
	UNUSED(k1);

#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
#endif
	unsigned int n = dat->n_observations;
	unsigned int p = dat->n_coordinates;
	data_t **nseeds = malloc(K * sizeof *nseeds);
	data_t *save_seed = malloc(p * sizeof *save_seed);
	unsigned int *nsd_idx = malloc(K * sizeof *nsd_idx);
	unsigned int save_sd_idx = 0;
	unsigned int k = 0;
	double hdis[K], dsum = 0;

	if (!nseeds || !nsd_idx || !save_seed)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nseeds");

	/* compute hamming distance of each cluster */
	for (k = 0; k < K; ++k)
		hdis[k] = 0;
	
	// if ALGORITHM FOUND A NULL CLUSTER
	int count = 0;
	for (k = 0; k < K; ++k)
		for (int j = 0; j < dat->n_coordinates; ++j)
			if (dat->best_modes[k][j] == 0)
				count++;
	
	for (unsigned int i = 0; i < n; ++i)
		if (count == K * dat->n_coordinates) {
			mmessage(WARNING_MSG, INTERNAL_ERROR, "null cluster?");
			hdis[dat->cluster_id[i]] += hd(dat->dmat[i],
					dat->seeds[dat->cluster_id[i]], p);
		} else
			hdis[dat->best_cluster_id[i]] += hd(dat->dmat[i],
				dat->best_modes[dat->best_cluster_id[i]], p);
	
	int null_cl = -1;
	/* choose seed to resample */
	for (k = 0; k < K; ++k) {
		if (count == K * dat->n_coordinates) {
			memcpy(seeds[k], dat->seeds[k], p * sizeof(**seeds));
			hdis[k] /= dat->cluster_size[k];
			if (dat->cluster_size[k] == 0)
				null_cl = k;
		} else {
			memcpy(seeds[k], dat->best_modes[k], p * sizeof(**seeds));
			hdis[k] /= dat->best_cluster_size[k];
		}
		
		dsum += hdis[k];
		nseeds[k] = seeds[k];
		nsd_idx[k] = sd_idx[k];
	}
	
	double r = runif(0, dsum);

	for (k = 0, dsum = 0; k < K && r > dsum; dsum += hdis[k++]);
	if (k)
		k--;
	if (null_cl != -1)
		k = null_cl;
	
	if (k != K - 1) {
		nseeds[K-1] = seeds[k];
		nseeds[k] = seeds[K-1];	/* we'll keep this one */
		nsd_idx[K-1] = sd_idx[k];
		nsd_idx[k] = sd_idx[K-1];	/* this one */
		save_sd_idx = sd_idx[K-1];
	}
	memcpy(save_seed, seeds[K-1], p*sizeof(*save_seed));

	do {
		initialize_proportional_to_abundance(dat, K, K-1, nseeds, nsd_idx);
	} while (!hd(nseeds[K-1], save_seed, p));	/* insist on new seed */

#ifdef DEBUGGING_CODE
	debug_msg(DEBUG_I >= fxn_debug, DEBUG_I, "Replacing %uth seed index %u "
		"with seed index %u\n", k, sd_idx[k], nsd_idx[K-1]);
#endif

	if (k != K - 1) {
		/* place new seed at old location */
		seeds[k] = nseeds[K-1];
		sd_idx[k] = nsd_idx[K-1];
		/* restore saved seed */
		memcpy(seeds[K-1], save_seed, p*sizeof(*save_seed));
		sd_idx[K-1] = save_sd_idx;
	}
   
	free(nseeds);
	free(nsd_idx);
	free(save_seed);

	return NO_ERROR;
} /* perturb_by_hd */


/**
 *
 * @param dat		data object (k-modes)
 * @param K		number of clusters
 * @param k1		unused
 * @param seeds		the current seeds
 * @param sd_idx	the current seed indices
 * @return		error status
 */
static int
perturb_methods(
				data *dat,
				unsigned int K,
				unsigned int k1,
				data_t **seeds,
				unsigned int *sd_idx,
				options *opt)
{																		   /**/
	UNUSED(k1);
	int err = NO_ERROR;
#ifdef DEBUGGING_CODE
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
#endif
	unsigned int p = dat->n_coordinates;
	data_t **nseeds = malloc(K * sizeof *nseeds);
	data_t *save_seed = malloc(p * sizeof *save_seed);
	unsigned int *nsd_idx = malloc(K * sizeof *nsd_idx);
	unsigned int save_sd_idx = 0;
	unsigned int k = 0;

	if (!nseeds || !nsd_idx || !save_seed)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nseeds");
	
	if (opt->perturb_selection == KMODES_PERTURB_HD) {
		k = select_by_hd (dat, seeds, nseeds, nsd_idx, sd_idx, K);
	}
	
	if (k != K - 1) {
		nseeds[K-1] = seeds[k];
		nseeds[k] = seeds[K-1];	/* we'll keep this one */
		nsd_idx[K-1] = sd_idx[k];
		nsd_idx[k] = sd_idx[K-1];	/* this one */
		save_sd_idx = sd_idx[K-1];
	}
	memcpy(save_seed, seeds[K-1], p*sizeof(*save_seed));

	do {
		if (opt->init_method == KMODES_INIT_AV07) {
			err = kmodes_init_av07(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, nseeds, nsd_idx, 0);
		} else if (opt->init_method == KMODES_INIT_AV07_GREEDY) {
			err = kmodes_init_av07(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, nseeds, nsd_idx, 1);
		} else if (opt->init_method == KMODES_INIT_CLB09) {
			err = kmodes_init_clb09(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, 0, nseeds, nsd_idx, UNSIGNED_INT);
		} else if (opt->init_method == KMODES_INIT_CLB09_RANDOM) {
			err = kmodes_init_clb09(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, 1, nseeds, nsd_idx, UNSIGNED_INT);
		} else if (opt->init_method == KMODES_INIT_H97) {
			err = kmodes_init_h97(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, 0, nseeds, nsd_idx);
		} else if (opt->init_method == KMODES_INIT_H97_RANDOM) {
			err = kmodes_init_h97(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, 1, nseeds, nsd_idx);
		} else if (opt->init_method == KMODES_INIT_HD17) {
			err = kmodes_init_hd17(dat->dmat, dat->n_observations, dat->n_coordinates, K, K-1, 0, nseeds, nsd_idx);
		} else
			replace_by_random(dat, K, K-1, nseeds, nsd_idx);
	} while (!hd(nseeds[K-1], save_seed, p));	/* insist on new seed */

#ifdef DEBUGGING_CODE
	debug_msg(DEBUG_I >= fxn_debug, DEBUG_I, "Replacing %uth seed index %u "
		"with seed index %u\n", k, sd_idx[k], nsd_idx[K-1]);
#endif

	if (k != K - 1) {
		/* place new seed at old location */
		seeds[k] = nseeds[K-1];
		sd_idx[k] = nsd_idx[K-1];
		/* restore saved seed */
		memcpy(seeds[K-1], save_seed, p*sizeof(*save_seed));
		sd_idx[K-1] = save_sd_idx;
	}
   
	free(nseeds);
	free(nsd_idx);
	free(save_seed);

	return err;
} /* perturb_methods */


int select_by_hd(data *dat, data_t **seeds,
				   data_t **nseeds, unsigned int *nsd_idx,
				   unsigned int *sd_idx, unsigned int K)
{
	unsigned int n = dat->n_observations;
	unsigned int p = dat->n_coordinates;
	double hdis[K], dsum = 0;
	unsigned int k = 0;
	
	/* compute hamming distance of each cluster */
	for (k = 0; k < K; ++k)
		hdis[k] = 0;
	
	// if algorithm gets a NULL cluster
	int count = 0;
	for (k = 0; k < K; ++k)
		for (int j = 0; j < dat->n_coordinates; ++j)
			if (dat->best_modes[k][j] == 0)
				count++;
	
	for (unsigned int i = 0; i < n; ++i)
		if (count == K * dat->n_coordinates) {
			mmessage(WARNING_MSG, INTERNAL_ERROR, "null cluster?");
			hdis[dat->cluster_id[i]] += hd(dat->dmat[i],
					dat->seeds[dat->cluster_id[i]], p);
		} else
			hdis[dat->best_cluster_id[i]] += hd(dat->dmat[i],
				dat->best_modes[dat->best_cluster_id[i]], p);
	
	int null_cl = -1;
	/* choose seed to resample */
	for (k = 0; k < K; ++k) {
		if (count == K * dat->n_coordinates) {
			memcpy(seeds[k], dat->seeds[k], p * sizeof(**seeds));
			hdis[k] /= dat->cluster_size[k];
			if (dat->cluster_size[k] == 0)
				null_cl = k;
		} else {
			memcpy(seeds[k], dat->best_modes[k], p * sizeof(**seeds));
			hdis[k] /= dat->best_cluster_size[k];
		}
		
		dsum += hdis[k];
		nseeds[k] = seeds[k];
		nsd_idx[k] = sd_idx[k];
	}
	
	double r = runif(0, dsum);

	for (k = 0, dsum = 0; k < K && r > dsum; dsum += hdis[k++]);
	if (k)
		k--;
	if (null_cl != -1)
		k = null_cl;
	
	return k;
} /* select_by_hd */

void replace_by_random(data *dat,
				  unsigned int K,
				  unsigned int k1,
				  data_t **seeds,
				  unsigned int *sd_idx) {
	data_t **x = dat->dmat;
	unsigned int p = dat->n_coordinates;
	unsigned int n = dat->n_observations;
	int repeat = 0;
	
	for (unsigned int k = k1; k < K; ++k) {
		if (repeat)
			--k;
		repeat = 0;
		
		size_t s = (unsigned int) runif(0, n);
		// get the center of obs with the chosen masked obs
	   
		sd_idx[k] = s;

		for (unsigned int l = 0; l < k1; ++l)
			if (!hd(seeds[l], x[s], p)) {
				repeat = 1;
				break;
			}
		memcpy(seeds[k], x[s], p * sizeof(**seeds));
	}
}

/**
 * Perturb for by masked abundance.
 *
 * @param dat	data object
 * @param opt	options object
 * @return	error status
 */
int perturb_general(data *dat, options *opt)
{
	return perturb_methods(dat, opt->K,
		opt->n_seed_set, dat->seeds, dat->seed_idx, opt);
} /* perturb */
