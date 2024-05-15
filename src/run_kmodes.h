/**
 * @file run_kmodes.h
 * @author Karin S. Dorman
 *
 * Header file for kmodes structs, extern functions and defines.
 */

#ifndef __H_RUN_KMODES__
#define __H_RUN_KMODES__

typedef struct _options options;
typedef struct _data data;

#include "kmodes.h"
#include "hash.h"

/**
 * Run options.
 */
struct _options {
	/* model */
	unsigned int K;		/*<! number of clusters */

	/* data */
	char const *datafile;		/*<! name of datafile */
	char const *data_outfile;	/*<! name of datafile to write */
	int subtract_one;		/*<! subtract one (to make 0-base) */
	unsigned int true_column;	/*<! supervised data: column of truth */
	unsigned int *true_cluster;	/*<! true cluster assignments */
	unsigned int *true_cluster_size;/*<! true cluster sizes */
	unsigned int true_K;		/*<! true number of clusters */
	data_t **true_modes;		/*<! true modes */
	unsigned int min_k;		/*<! minimum K */
	unsigned int max_k;		/*<! maximum K */
	unsigned int n_k;		/*<! number of K [min_k, max_k] */
	unsigned int min_nk;		/*<! minimum K for no update ini */
	unsigned int max_nk;		/*<! minimum K for no update ini */
	unsigned int n_nk;		/*<! number of no up K */

	/* run conditions */
	kmodes_options *kopt;	/*<! kmodes run options passed to algorithms */
	unsigned int n_max_iter;/*<! max. number of iterations */
	uint_fast64_t seed;	/*<! random number seed */
	uint_fast64_t seed2;	/*<! second random number seed (set_seed() needs 2) */
	int kmodes_algorithm;	/*<! run Lloyd's, Huang's or HW k-modes */
	int use_hartigan;	/*<! use Hartigan updates (Lloyd's | Huang's) */
	int update_modes;	/*<! update modes */
	int use_qtran;		/*<! use quick-transfer phase */
	int weight;		/*<! weighting method (not implemented) */
	double seconds;		/*<! number of seconds to initialize */
	double target;		/*<! target minimum: best from prev. run */
	int shuffle;		/*<! shuffle input order */
	int continue_run;	/*<! continue previous run */
	int select_k;		/*<! select k */
	int dm_method;		/*<! select K by DM method */
	double n_effective_coordinates;	/*<! effective no. indep. coordinates */
	char const ***result_files;	/*<! names of output files per K */
	unsigned int *n_result_files;	/*<! number of output files per K */
	char const ***nup_ini_files;	/*<! ini file names w/o update per K */
	unsigned int *n_nup_ini_files;	/*<! no. ini files w/o update per K */
	unsigned int n_bootstrap;	/*<! no. of bootstraps in K-select */

	/* initialization */
	unsigned int n_init;	/*<! number of random initializations */
	unsigned int n_inner_init;	/*<! number of "inner" initializations */
	int init_method;	/*<! initialization method to use */
	unsigned int *seed_idx;	/*<! user-provided seed indices */
	char const *pfile;	/*<! partition file */
	char const *sfile;	/*<! seed file */
	char const *mfile;	/*<! mode file */
	char const *mfile_out;	/*<! name of mode file to output */
	unsigned int n_sd_idx;	/*<! tmp: length of seed_idx */
	data_t **seed_set;	/*<! available seeds */
	unsigned int n_seed_set;/*<! number of seeds in seed set */
	int abunk;		/*<! number of masked abundance>1 */
	unsigned int resample_sites;	/*<! resample sites every int inits */
    int perturb_selection; /*<! perturbation selection method */
    unsigned int n_perturb;       /*<! perturb this many seeds when re-initialize */
    unsigned int inner_perturb; /*<! 2 settings for no. of perturbations per initialization */

	/* output */
	char const *ini_file;	/*<! initialization data outfile */
	char const *soln_file;	/*<! solution file */
	int info;		/*<! level of information to output */
	int quiet;		/*<! be quiet */

	/* simulation */
	int simulate;			/*<! request to simulate data */
	unsigned int sim_K;		/*<! number of simulated clusters */
					/* options:true_K <= options:sim_K */
	int require_sim_K;		/*<! require true_K == sim_K */
	double *sim_alpha;		/*<! dirichlet parameters [optional] */
	double *sim_pi;			/*<! mixing proportions */
	double sim_between_t;		/*<! between variation */
	double sim_within_t;		/*<! within variation */
	double sim_within_prob, sim_between_prob;
	unsigned int sim_n_observations;/*<! number of observations */
	unsigned int sim_n_coordinates;	/*<! number of coordinates */
	data_t sim_n_categories;	/*<! number of categories */
	data_t **sim_modes;		/*<! simulation (or true) modes */
	unsigned int *sim_cluster;	/*<! simulated cluster assignments */
}; /* options */


/**
 * Data.  Store the input data.
 */
struct _data {
	data_t *data;			/*<! data */
	data_t **dmat;			/*<! data as matrix */
	data_t **masked_dmat;		/*<! masked copy of data */
	unsigned int n_observations;	/*<! number of observations */
	unsigned int n_coordinates;	/*<! number of coordinates (all categorical) */
	data_t *n_categories;		/*<! number of categories in each coordinate */
	data_t max_n_categories;	/*<! max. number categories per coordinate */
	size_t tot_n_categories;	/*<! total number of categories */

	/* initialization */
	data_t **seeds;			/*<! chosen seeds */
	data_t **ini_seeds;		/*<! trial seeds */
	unsigned int *seed_idx;		/*<! seed indices */
	unsigned int *ini_seed_idx;	/*<! trial seed indices */
	unsigned int n_init;		/*<! number of initializations done */
	uint8_t use_ini;		/*<! using ini_* versions or not */
	unsigned int n_masked;		/*<! number of masked sites */
    
	/* hash table */
	hash *seq_count;		/*<! frequency of unique sequences (hash table) */
	unsigned int hash_length;	/*<! num of unique sequences in hash table */
	double *entropy;		/*<! per-coordinate entropy */
	double total_entropy;		/*<! sum of entropy */
	size_t *entropy_order;		/*<! coordinates sorted by entropy */
	unsigned int *mask;		/*<! for masked abundance */

	/* current solution */
	double total;			/*<! current criterion */
	unsigned int *cluster_id;	/*<! cluster assignments */
	unsigned int *obsn_idx;		/*<! used if shuffling */
	double *criterion;		/*<! criterion */
	unsigned int *cluster_size;	/*<! cluster sizes */
	unsigned int iter;		/*<! iterations */

	/* best solution */
	double best_total;		/*<! best total criterion */
	double best_rand;		/*<! if true clusters known, not */
					/*   necessarily at best objective */
	double best_rand_at_opt;	/*<! rand at first optimum */
	double best_mi_at_opt;		/*<! mi at first optimum */
	double best_vi_at_opt;		/*<! vi at first optimum */
	unsigned int best_total_cnt;	/*<! in case, same optimum repeated */
	unsigned int *best_seed_idx;	/*<! seeding at optimum */
	data_t **best_modes;		/*<! estimated modes at optimum */
	double *best_criterion;		/*<! criterion at optimum */
	unsigned int *best_cluster_id;	/*<! cluster assignments at optimum */
	unsigned int *best_cluster_size;/*<! cluster sizes at optimum */
	unsigned int *best_obsn_idx;	/*<! used if shuffling at optimum */

	/* summary statistics */
	double seconds;			/*<! seconds used */
	double avg_cost, sd_cost;	/*<! minimum criterion */
	double avg_iter, sd_iter;	/*<! iterations to convergence */
	double avg_ar, sd_ar;		/*<! adjusted rand (truth known) */
	double avg_mi, sd_mi;		/*<! mutual information (truth known) */
	double avg_vi, sd_vi;		/*<! variance of info (truth known) */
	double avg_time, sd_time;	/*<! inits to target */
	unsigned int ntimes;		/*<! times hit target */

	/* internal use */
	double uncounted_seconds;
	double first_cost;
	double worst_cost;
	unsigned int max_iter;
	unsigned int ctime;
}; /* data */

int make_options(options **opt);
int make_data(data **data, options *opt);
int allocate_data(data *dat, unsigned int no_data);
int finish_make_data(data *dat, options *opt);
int drop_empty_clusters(data *dat, options *opt);
int fix_categories(data *dat, options *opt);
int run_kmodes(data *dat, options *opt);
void free_options(options *opt);
void free_data(data *dat);


#endif
