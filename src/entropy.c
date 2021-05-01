
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>

#include "run_kmodes.h"
#include "array.h"
#include "error.h"
#include "io.h"
#include "math.h"
#include "entropy.h"
#include "order.h"
#include "hash.h"

static inline double hd(data_t *x, data_t *y, unsigned int p);
static int perturb_by_hd(data *dat, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
static int initialize_proportional_to_abundance(data *dat, unsigned int K, unsigned int k1, data_t **seeds,
    unsigned int *sd_idx);
static double entropy(data_t *x, int unique_c, int len);
static double * entropy_calc(data_t **x, unsigned int n, unsigned int p, data_t max_n_categories);
static unsigned int mask_idx(double *entr, unsigned int p, double persentage, size_t **indices);
int binarysearch(data_t **copied, data *dat, size_t *indices,
                 int l, int r, unsigned int K, int abunk);

static inline double hd(data_t *x, data_t *y, unsigned int p)
{
    double d = 0;
    for (unsigned int j = 0; j < p; ++j)
        d += (x[j] != y[j]);
    return d;
} /* hd */

static double entropy(data_t *x, int unique_c, int len) {
    double entro = 0;
    int hist[unique_c];
    /* assume categories are 0, 1, 2, ... without skips */
    for(int i = 0; i < unique_c; ++i)
        hist[i] = 0;
    for(int i = 0; i < len; ++i)
        hist[x[i]]++;
    
    for(int i = 0; i < unique_c; ++i) {
        if (hist[i] == 0)
            continue;
        entro -= ((double)hist[i] / len) * (log2((double)hist[i]) - log2(len));
    }
    
    return entro;
}

static double * entropy_calc(data_t **x, unsigned int n, unsigned int p, data_t max_n_categories) {
    unsigned int i, j;
    double *entr = NULL;
    entr = calloc(sizeof(*entr), p);
    fprintf(stderr, "Entropy: ");
    for (j = 0; j < p; ++j) {
        data_t tmp[n];
        for (i = 0; i < n; ++i)
            tmp[i] = x[i][j];
        entr[j] = entropy(tmp, max_n_categories, n);
        fprintf(stderr, "%lf ", entr[j]);
    }
    fprintf(stderr, "\n");
    return entr;
}

static unsigned int mask_idx(double *entr, unsigned int p, double persentage, size_t **indices) {
//    size_t *pos = NULL;
    UNUSED(persentage);
//    *indices = malloc(sizeof(size_t) * p);
    
    *indices = qorder_double(entr, p); // get the order of entr
    
    for (unsigned int j = 0, i = p - 1; j < i; ++j, --i) {
        size_t tmp = (*indices)[j];
        (*indices)[j] = (*indices)[i];
        (*indices)[i] = tmp;
    }
    fprintf(stderr, "Ranked sites: ");
    for (unsigned int j = 0; j < p; ++j)
        fprintf(stderr, "%zu ", (*indices)[j]);
    fprintf(stderr, "\n");
//    free(pos);
    return NO_ERROR;
}

int binarysearch(data_t **copied, data *dat, size_t *indices,
                 int l, int r, unsigned int K, int abunk) {
    unsigned int i, j, m;
    unsigned int mid = 0;
    if(r > l) {
        mid = l + (r - l) / 2;
        // this is not efficient
        for (j = 0; j < mid; ++j)
            for (i = 0; i < dat->n_observations; ++i)
                copied[i][indices[j]] = 0;
        
        /* build hash table */
        dat->hash_length = 0;
        dat->seq_count = NULL;
        for (i = 0; i < dat->n_observations; ++i)
            dat->hash_length += add_sequence(&dat->seq_count, copied[i],
                                             dat->n_coordinates, i);
        m = 0;
        int pool = 0;
        for (hash *s = dat->seq_count; s != NULL; s = s->hh.next, ++m) {
            if(s->count > 1)
                pool++;
        }
        if(pool <= abunk)
            return mid;
//        else if(pool < abunk) {// mask less, revert copied first
//            for (j = 0; j < mid; ++j)
//                for (i = 0; i < dat->n_observations; ++i)
//                    copied[i][indices[j]] = dat->dmat[i][indices[j]];
//            return binarysearch(copied, dat, indices, l, mid - 1, K, abunk);
//        }
        else
            return binarysearch(copied, dat, indices, mid + 1, r, K, abunk);
    } else {
        mmessage(WARNING_MSG, CUSTOM_ERROR, "Please give a large -a option.\n");
        dat->hash_length = 0;
        dat->seq_count = NULL;
        for (i = 0; i < dat->n_observations; ++i)
            dat->hash_length += add_sequence(&dat->seq_count, dat->dmat[i],
                                             dat->n_coordinates, i);
        return mid;
    }
}

int mask_nhash(data *dat, unsigned int K, int abunk) {
    int err = NO_ERROR;
    unsigned int p = dat->n_coordinates, n = dat->n_observations;
    unsigned int i, j, m;
    data_t **x = dat->dmat;
    size_t *indices = NULL;
    
    double *entropy = entropy_calc(x, n, p , dat->max_n_categories);
    err = mask_idx(entropy, p, 1, &indices);
    
    data_t **copied = NULL;
    copied = malloc(n * sizeof *x);
    for (i = 0; i < n; ++i)
        copied[i] = malloc(p * sizeof **x);
    
    for (i = 0; i < n; ++i)
        for (j = 0; j < p; ++j)
            copied[i][j] = x[i][j];
    // use binary search to find out the number of masked position
    unsigned int topk = binarysearch(copied, dat, indices, 0, p - 1, K, abunk);
    fprintf(stderr, "Masked sites: ");
    for (j = 0; j < topk; ++j)
        fprintf(stderr, "%zu ", indices[j]);
    fprintf(stderr, "\n");

    m = 0;
//    for (hash *s = dat->seq_count; s != NULL; s = s->hh.next, ++m)
//        fprintf(stderr, "Abundance[%zu(%u)=", s->idx, m);
    /* store index of reads for all unique sequences */
    for (i = 0; i < n; ++i)
        if ((err = add_seq_idx(dat->seq_count, copied[i], p, i)))
            return err;

    for (hash *s = dat->seq_count; s != NULL; s = s->hh.next, ++m) {
        fprintf(stderr, "Abundance[%zu(%u)=", s->idx, m);
        for (unsigned int j = 0; j < p; ++j)
            fprintf(stderr, "%d", copied[s->idx][j]);
        fprintf(stderr, "]: %u\n", s->count);
    }
    
    free(indices);
    free(copied);
    return err;
}

/* Init with abundance by entropy mask high entropy site */
int
initialize_by_abundance(
                data *dat,
                unsigned int K,
                unsigned int k1,
                data_t **seeds,
                unsigned int *sd_idx)
{
    int err = NO_ERROR;
    unsigned int p = dat->n_coordinates;
    unsigned int j, k, l, nties = 0;
    data_t **x = dat->dmat;
    hash *sc = dat->seq_count, *snext = NULL;
    
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
                    memcpy(seeds[j], x[sc->idx], p * sizeof(**seeds));
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
            k = (unsigned int) nties * (rand() / (RAND_MAX + 1.));

            if (used[k])
                continue;

            /* forward to it */
            snext = sc;
            for (l = 0; l < k; ++l)
                snext = snext->hh.next;

            memcpy(seeds[j], x[snext->idx], p * sizeof(**seeds));
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

int perturbs(data *dat, options *opt)
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
        r = dsum * rand() / (RAND_MAX + 1.);

        for (useq = dat->seq_count, dsum = useq->count; useq != NULL && r > dsum; useq = useq->hh.next, dsum += useq ? useq->count : 0);

        s = useq->idx;
        // get the center of obs with the chosen masked obs
        memcpy(seeds[k], x[s], p * sizeof(**seeds));
        sd_idx[k] = s;

        for (unsigned int l = 0; l < k1; ++l)
            if (!hd(seeds[l], x[s], p)) {
                repeat = 1;
                break;
            }
        
//        s = dat->cluster_id[s];
//        memcpy(seeds[k], dat->seeds[s], p * sizeof(**seeds));
//        sd_idx[k] = s;
//
//        for (unsigned int l = 0; l < k1; ++l)
//            if (!hd(seeds[l], dat->seeds[s], p)) {
//                repeat = 1;
//                break;
//            }
    }

    return NO_ERROR;
} /* initialize_proportional_to_abundance */

/**
 *
 * @param dat        data object (k-modes)
 * @param K        number of clusters
 * @param k1        unused
 * @param seeds        the current seeds
 * @param sd_idx    the current seed indices
 * @return        error status
 */
static int
perturb_by_hd(
                data *dat,
                unsigned int K,
                unsigned int k1,
                data_t **seeds,
                unsigned int *sd_idx)
{
    UNUSED(k1);

    unsigned int n = dat->n_observations;
    unsigned int p = dat->n_coordinates;
    data_t **nseeds = malloc(K * sizeof *nseeds);
    data_t *save_seed = malloc(p * sizeof *save_seed);
    unsigned int *nsd_idx = malloc(K * sizeof *nsd_idx);
    unsigned int save_sd_idx = 0;
    unsigned int k = 0;
    double ll[K], dsum = 0;

    if (!nseeds || !nsd_idx)
        return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "nseeds");

    /* compute hamming distance of each cluster */
    for (k = 0; k < K; ++k)
        ll[k] = 0;
    
    // if ALGORITHM FOUND A NULL CLUSTER
    int count = 0;
    for (k = 0; k < K; ++k)
        for(int j = 0; j < dat->n_coordinates; ++j)
            if(dat->best_modes[k][j] == 0)
                count++;
    
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < p; ++j) {
            if (count == K * dat->n_coordinates) {
                if (dat->dmat[i][j]
                    != dat->seeds[dat->cluster_id[i]][j])
                    ll[dat->cluster_id[i]]
                            += 1;
            } else {
                if (dat->dmat[i][j]
                    != dat->best_modes[dat->best_cluster_id[i]][j])
                    ll[dat->best_cluster_id[i]]
                            += 1; }
        }
    }
    
    int null_cl = -1;
    /* choose seed to resample */
    for (k = 0; k < K; ++k) {
        if (count == K * dat->n_coordinates) {
            memcpy(seeds[k], dat->seeds[k], p * sizeof(**seeds));
            ll[k] /= dat->cluster_size[k];
            if(dat->cluster_size[k] == 0)
                null_cl = k;
        } else {
            memcpy(seeds[k], dat->best_modes[k], p * sizeof(**seeds));
            ll[k] /= dat->best_cluster_size[k];
        }
        
        dsum += ll[k];
        nseeds[k] = seeds[k];
        nsd_idx[k] = sd_idx[k];
    }
    
    double r = dsum * rand() / (RAND_MAX + 1.);

    for (k = 0, dsum = 0; k < K && r > dsum; dsum += ll[k++]);
    if (k)
        k--;
    if(null_cl != -1)
        k = null_cl;
    
    if (k != K - 1) {
        nseeds[K-1] = seeds[k];
        nseeds[k] = seeds[K-1];    /* we'll keep this one */
        nsd_idx[K-1] = sd_idx[k];
        nsd_idx[k] = sd_idx[K-1];    /* this one */
        save_sd_idx = sd_idx[K-1];
    }
    memcpy(save_seed, seeds[K-1], p*sizeof(*save_seed));

    do {
        initialize_proportional_to_abundance(dat, K, K-1, nseeds, nsd_idx);
    } while (!hd(nseeds[K-1], save_seed, p));    /* insist on new seed */

    //fprintf(stderr, "Choosing new seed %u: %u\n", k, sd_idx[k]);

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
