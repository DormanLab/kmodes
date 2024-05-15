/**
 * @file run_kmodes_r.c
 *
 * Purpose: C wrapper for k-modes.  R wrapper in run_kmodes.R.
 */

#include <R.h>
#include <Rinternals.h>

#include "run_kmodes.h"
#include "error.h"

enum RETURN_SLOT {
	TOTAL_CRITERION_SLOT,
	CRITERIA_SLOT,
	CLUSTER_SIZE_SLOT,
	PARTITION_SLOT,
	MODES_SLOT,
	AVERAGE_CRITERION_SLOT,
	NUMBER_INITIALIZATIONS_SLOT,
	BEST_RAND_SLOT,
	MAX_RAND_SLOT,
	AVERAGE_RAND_SLOT,
	NUMBER_SLOTS,
};

char const *slot_names[NUMBER_SLOTS] = {
	"criterion", "cluster.criteria",
	"cluster.sizes", "partition",
	"modes", "average.criterion",
	"number.initializations",
	"best.ari", "max.ari", "average.ari"
};

int load_data(data *dat, options *opt, SEXP data_r, int *);

/**
 * This function is called from R.
 *
 * @param data_r	matrix of integer-interpretable types (pxn INTSXP -- transposed in R)
 * @param K_r		number of clusters (INTSXP scalar)
 * @param algorithm_r	algorithm (STRSXP scalar)
 * @param init_method_r	initialization method (STRSXP scalar)
 * @param ninit_r	number of initializations (INTSXP scalar)
 * @param true_clus_r	true cluster (n-length INTSXP vector)
 * @param shuffle_r	shuffle the observation order in data
 * @param verbosity_r	verbosity (0, 1, 2, 3)
 * @return		list
 */
SEXP run_kmodes_r(	SEXP data_r,
			SEXP K_r,
			SEXP algorithm_r,
			SEXP init_method_r,
			SEXP ninit_r,
			SEXP true_clus_r,
			SEXP shuffle_r,
			SEXP verbosity_r)
{
	int err = NO_ERROR;
	int nprotect = 0;
	data *dat = NULL;
	options *opt = NULL;
	int *iptr;
	SEXP return_list = R_NilValue;
	SEXP best_modes;
	SEXP return_list_names;

	if ((err = make_options(&opt)))
		goto EXIT_RUN_KMODES_R;

	if (!isInteger(K_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type "
						"for number of clusters.\n");
		goto EXIT_RUN_KMODES_R;
	}
	if (!isInteger(ninit_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type "
					"for number of initializations.\n");
		goto EXIT_RUN_KMODES_R;
	}
	if (!isString(algorithm_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type "
				"for k-modes algorithm, option algorithm.\n");
		goto EXIT_RUN_KMODES_R;
	}
	if (!isString(init_method_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type "
			"for initialization method, option init.method.\n");
		goto EXIT_RUN_KMODES_R;
	}
	if (!isLogical(shuffle_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type for shuffle.\n");
		goto EXIT_RUN_KMODES_R;
	}
	if (!isInteger(verbosity_r)) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Invalid type for verbosity level.\n");
		goto EXIT_RUN_KMODES_R;
	}

	/* copy user-supplied options into options object;
	 * not much is available yet
	 */

	opt->K = * INTEGER(K_r);
	opt->n_init = * INTEGER(ninit_r);
	opt->quiet = * INTEGER(verbosity_r);
	opt->shuffle = * LOGICAL(shuffle_r);

	char const *alg = CHAR(STRING_ELT(algorithm_r, 0));
	char const *ini = CHAR(STRING_ELT(init_method_r, 0));

	if (!strcmp(alg, "h97")) {
		opt->kmodes_algorithm = KMODES_HUANG;
	} else if (!strcmp(alg, "lloyd")) {
		opt->kmodes_algorithm = KMODES_LLOYD;
	} else if (!strcmp(alg, "hw")) {
		opt->kmodes_algorithm = KMODES_HARTIGAN_WONG;
	} else {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Unrecognized "
			"algorithm '%s', option algorithm.\n", alg);
		goto EXIT_RUN_KMODES_R;
	}

	if (!strcmp(ini, "rnd")) {
		opt->init_method = KMODES_INIT_RANDOM_SEEDS;
	} else if (!strcmp(ini, "h97")) {
		opt->init_method = KMODES_INIT_H97;
	} else if (!strcmp(ini, "h97rnd")) {
		opt->init_method = KMODES_INIT_H97_RANDOM;
	} else if (!strcmp(ini, "hd17")) {
		opt->init_method = KMODES_INIT_HD17;
	} else if (!strcmp(ini, "clb09")) {
		opt->init_method = KMODES_INIT_CLB09;
	} else if (!strcmp(ini, "clb09rnd")) {
		opt->init_method = KMODES_INIT_CLB09_RANDOM;
	} else if (!strcmp(ini, "av07")) {
		opt->init_method = KMODES_INIT_AV07;
	} else if (!strcmp(ini, "av07grd")) {
		opt->init_method = KMODES_INIT_AV07_GREEDY;
	} else if (!strcmp(ini, "rndp")) {
		opt->init_method = KMODES_INIT_RANDOM_FROM_PARTITION;
	} else if (!strcmp(ini, "rnds")) {
		opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
	/* NOT YET: KMODES_INIT_USER_SEEDS */
	} else {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Unrecognized "
					"initialization method '%s'.\n", ini);
		goto EXIT_RUN_KMODES_R;
	}

	if (!isNull(true_clus_r)) {
		if (!isInteger(true_clus_r)) {
			err = INVALID_USER_INPUT;
			mmessage(ERROR_MSG, err, "Invalid type for true "
						"cluster membership.\n");
			goto EXIT_RUN_KMODES_R;
		}

		opt->true_cluster = (unsigned int *)INTEGER(true_clus_r);

		debug_msg(opt->quiet >= VERBOSE, opt->quiet, "True cluster mem"
			"bership for %d observations.\n", length(true_clus_r));
	}

	/* create data object */
	if ((err = make_data(&dat, opt))) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Failed to allocate data object.\n");
		goto EXIT_RUN_KMODES_R;
	}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Made data object.\n");

	/* load data into data object */
	if ((err = load_data(dat, opt, data_r, &nprotect))) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Failed to load data from R.\n");
		goto EXIT_RUN_KMODES_R;
	}

	debug_msg(opt->quiet >= DEBUG_I, opt->quiet, "Loaded data object of "
		"dimension %dx%d.\n", dat->n_observations, dat->n_coordinates);

	/* shuffle overwrites options::true_cluster */
	if (opt->true_cluster && opt->shuffle) {
		opt->true_cluster = malloc(dat->n_observations
						* sizeof(*opt->true_cluster));
		if (!opt->true_cluster) {
			err = MEMORY_ALLOCATION;
			mmessage(ERROR_MSG, err, "options::true_cluster\n");
			goto EXIT_RUN_KMODES_R;
		}
		int *itmp = INTEGER(true_clus_r);
		for (unsigned int i = 0; i < dat->n_observations; ++i)
			opt->true_cluster[i] = itmp[i];

/*
		memcpy(opt->true_cluster, (unsigned int *)INTEGER(true_clus_r),
				dat->n_observations * sizeof(*opt->true_cluster));
*/
	}

	/* finish preparing data object */
	if ((err = finish_make_data(dat, opt))) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "Failed to finalize data object.\n");
		goto EXIT_RUN_KMODES_R;
	}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Data finished.\n");

	/* run k-modes: information returned in dat */
	if ((err = run_kmodes(dat, opt))) {
		err = INVALID_USER_INPUT;
		mmessage(ERROR_MSG, err, "K-modes failed.\n");
		goto EXIT_RUN_KMODES_R;
	}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Completed k-modes runs.\n");

	/* set up return object */
	PROTECT(return_list = allocVector(VECSXP, NUMBER_SLOTS));
	++nprotect;

	if (!return_list) {
		err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, err, "return_list");
		goto EXIT_RUN_KMODES_R;
	}

	SET_VECTOR_ELT(return_list, TOTAL_CRITERION_SLOT,
					ScalarReal(dat->best_total));
	debug_msg(opt->quiet >= VERBOSE, opt->quiet, "Best criterion: %f\n",
								dat->best_total);

	SET_VECTOR_ELT(return_list, CRITERIA_SLOT,
					allocVector(REALSXP, opt->K));
	memcpy(REAL(VECTOR_ELT(return_list, CRITERIA_SLOT)),
		dat->best_criterion, opt->K * sizeof(*dat->best_criterion));
	debug_msg(opt->quiet >= DEBUG_I, opt->quiet, "Criteria of first three "
				"clusters: %f %f %f\n", dat->best_criterion[0],
				dat->best_criterion[1], dat->best_criterion[2]);

	SET_VECTOR_ELT(return_list, CLUSTER_SIZE_SLOT,
					allocVector(INTSXP, opt->K));
	memcpy(INTEGER(VECTOR_ELT(return_list, CLUSTER_SIZE_SLOT)),
					dat->best_cluster_size, opt->K
					* sizeof(*dat->best_cluster_size));
	debug_msg(opt->quiet >= DEBUG_I, opt->quiet, "First three cluster "
				"sizes: %d %d %d\n", dat->best_cluster_size[0],
			dat->best_cluster_size[1], dat->best_cluster_size[2]);

	SET_VECTOR_ELT(return_list, PARTITION_SLOT,
				allocVector(INTSXP, dat->n_observations));
	memcpy(INTEGER(VECTOR_ELT(return_list, PARTITION_SLOT)),
				dat->best_cluster_id, dat->n_observations
					* sizeof(*dat->best_cluster_id));
	debug_msg(opt->quiet >= DEBUG_III, opt->quiet,
						"Added best cluster IDs.\n");

	PROTECT(best_modes = allocMatrix(INTSXP, opt->K, dat->n_coordinates));
	++nprotect;
	iptr = INTEGER(best_modes);
	for (unsigned int k = 0; k < opt->K; ++k)
		for (unsigned int j = 0; j < dat->n_coordinates; ++j)
			iptr[opt->K*j + k] = dat->best_modes[k][j];
	SET_VECTOR_ELT(return_list, MODES_SLOT, best_modes);
	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Added best modes.\n");

	if (dat->n_init) {
		SET_VECTOR_ELT(return_list, AVERAGE_CRITERION_SLOT,
						ScalarReal(dat->avg_cost));
		SET_VECTOR_ELT(return_list, NUMBER_INITIALIZATIONS_SLOT,
						ScalarInteger(dat->n_init));
		if (opt->simulate || opt->true_cluster) {
			SET_VECTOR_ELT(return_list, BEST_RAND_SLOT,
						ScalarReal(dat->best_rand_at_opt));
			SET_VECTOR_ELT(return_list, MAX_RAND_SLOT,
						ScalarReal(dat->best_rand));
			SET_VECTOR_ELT(return_list, AVERAGE_RAND_SLOT,
						ScalarReal(dat->avg_ar));
		} else {
			SET_VECTOR_ELT(return_list, BEST_RAND_SLOT,
							ScalarReal(NA_REAL));
			SET_VECTOR_ELT(return_list, MAX_RAND_SLOT,
							ScalarReal(NA_REAL));
			SET_VECTOR_ELT(return_list, AVERAGE_RAND_SLOT,
							ScalarReal(NA_REAL));
		}
	} else {
		SET_VECTOR_ELT(return_list, AVERAGE_CRITERION_SLOT,
						ScalarReal(NA_REAL));
		SET_VECTOR_ELT(return_list, NUMBER_INITIALIZATIONS_SLOT,
						ScalarInteger(NA_INTEGER));
		SET_VECTOR_ELT(return_list, BEST_RAND_SLOT,
							ScalarReal(NA_REAL));
		SET_VECTOR_ELT(return_list, MAX_RAND_SLOT,
							ScalarReal(NA_REAL));
		SET_VECTOR_ELT(return_list, AVERAGE_RAND_SLOT,
							ScalarReal(NA_REAL));
	}
	debug_msg(opt->quiet >= DEBUG_III, opt->quiet,
					"Added average information.\n");

	/* give the list members names */
	PROTECT(return_list_names = allocVector(STRSXP, NUMBER_SLOTS));
	++nprotect;
	SET_STRING_ELT(return_list_names, TOTAL_CRITERION_SLOT,
				mkChar(slot_names[TOTAL_CRITERION_SLOT]));
	SET_STRING_ELT(return_list_names, CRITERIA_SLOT,
				mkChar(slot_names[CRITERIA_SLOT]));
	SET_STRING_ELT(return_list_names, CLUSTER_SIZE_SLOT,
				mkChar(slot_names[CLUSTER_SIZE_SLOT]));
	SET_STRING_ELT(return_list_names, PARTITION_SLOT,
					mkChar(slot_names[PARTITION_SLOT]));
	SET_STRING_ELT(return_list_names, MODES_SLOT,
					mkChar(slot_names[MODES_SLOT]));
	SET_STRING_ELT(return_list_names, AVERAGE_CRITERION_SLOT,
				mkChar(slot_names[AVERAGE_CRITERION_SLOT]));
	SET_STRING_ELT(return_list_names, NUMBER_INITIALIZATIONS_SLOT,
			mkChar(slot_names[NUMBER_INITIALIZATIONS_SLOT]));
	SET_STRING_ELT(return_list_names, BEST_RAND_SLOT,
					mkChar(slot_names[BEST_RAND_SLOT]));
	SET_STRING_ELT(return_list_names, AVERAGE_RAND_SLOT,
					mkChar(slot_names[AVERAGE_RAND_SLOT]));
	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Set strings.\n");
	setAttrib(return_list, R_NamesSymbol, return_list_names);

EXIT_RUN_KMODES_R:

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Freeing memory.\n");
	if (opt)
		free_options(opt);
	if (dat)
		free_data(dat);

	if (nprotect)
		UNPROTECT(nprotect);

	return return_list;
} /* run_kmodes_r */

/**
 * Load data from R into C data structure.
 *
 * @param dat		pointer to data object
 * @param opt		pointer to options object
 * @param data_r	matrix from R
 * @param nprotect	number of objects needing protection
 * @return		error status
 */
int load_data(data *dat, options *opt, SEXP data_r, int *nprotect)
{
	int err = NO_ERROR;
	int lnprotect = 0;
	int *iptr;
	unsigned int n, p;
	data_t *dptr;
	SEXP data_c = NULL;
	SEXP dim = getAttrib(data_r, R_DimSymbol);

	if (!isInteger(data_r)) {
		PROTECT(data_c = coerceVector(data_r, INTSXP));
		++lnprotect;
		iptr = INTEGER(data_c);
	} else {
		iptr = INTEGER(data_r);
	}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Coerced data.\n");

	/* recall data are transposed */
	p = dat->n_coordinates = * INTEGER(dim);
	n = dat->n_observations = INTEGER(dim)[1];

	debug_msg(opt->quiet >= VERBOSE, opt->quiet, "Read data of size "
							"%dx%d\n", n, p);

	if ((err = allocate_data(dat, 1)))
		goto EXIT_LOAD_DATA;

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Allocated data.\n");

	dptr = dat->data;

	/* convert R data to C data: column-major -> row-major */
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < p; ++j) {
			*dptr = iptr[i*p + j];
			if (opt->subtract_one)
				-- (*dptr);
			if (dat->n_categories[j] < *dptr + 1u)
				dat->n_categories[j] = *dptr + 1u;
			dptr++;
		}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Copied data.\n");

	/* recount and reindex clusters if any are empty */
	if (opt->true_cluster) {
		unsigned int min = 1;

		opt->true_K = 0;
		for (size_t i = 0; i < n; ++i) {
			if (opt->true_cluster[i] > opt->true_K)
				opt->true_K = opt->true_cluster[i];
			if (opt->true_cluster[i] < min)
				min = opt->true_cluster[i];
		}
		if (!min)	/* 0-based cluster indices */
			++opt->true_K;
		if ((err = drop_empty_clusters(dat, opt)))
			goto EXIT_LOAD_DATA;
	}

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet,
						"Dropped empty clusters.\n");

	/* fix assumption of 0, 1, 2, ... categories; reset category counts */
	if ((err = fix_categories(dat, opt)))
		goto EXIT_LOAD_DATA;

	debug_msg(opt->quiet >= DEBUG_III, opt->quiet, "Fixed categories.\n");
	*nprotect += lnprotect;

	return err;

EXIT_LOAD_DATA:

	if (lnprotect)
		UNPROTECT(lnprotect);

	return err;
} /* load_data */
