################################################################################
# @file run_kmodes.R
#
# Purpose: R wrapper for k-modes algorithms available in kmodes C software.
################################################################################

#' Categorical data clustering with k-modes
#'
#' Estimate clusters for categorical data using k-modes. Various k-modes
#' algorithms are available and several initialization methods.
#'
#' The k-modes objective function for hard clustering categorical data is the
#' sum across all K clusters of the Hamming distances between the cluster mode
#' and the cluster members.  The objective function can be optimized using many
#' algorithms, the most popular of which is Huang's algorithm.  Huang's algorithm
#' \cite{Huang1997b} is analogous to MacQueen's algorithm \cite{Macqueen1967} for
#' k-means.  Each observation is considered in turn.  If it is closer to the mode
#' of another cluster, it is moved to the new cluster and the modes of both
#' affected clusters are updated.  The process continues until there is no move
#' after one full cycle through the observations. Huang's algorithm is also
#' implemented in [klaR::kmodes()]. The implementation in this package [kmodes()]
#' is several times faster because the algorithm is implemented in C.
#' 
#' In addition to Huang's algorithm, there is Chaturvedi, Green, and Carroll's
#' algorithm \cite{Chaturvedi2001}, where modes are not updated until all
#' observations have been processed.  This algorithm is analogous to Lloyd's
#' algorithm \cite{Lloyd1982} for k-means.
#' 
#' This package also adapts Hartigan and Wong's efficient algorithm for k-means
#' \cite{Hartigan1979} to k-modes. In this algorithm, an observation is not moved
#' from one mode to another \emph{unless} the objective function will improve.
#' While more computationally expensive per iteration, the algorithm is known to
#' converge to the local optimum, whereas all previous algorithms may terminate
#' early. As a result, it usually obtains a better optimum after fewer
#' initializations and less time.
#' 
#' Because all these algorithms are local optimization algorithms, good
#' initialization or repeated random initialization is critical to find the global
#' optimum. There have been several initialization methods suggested for k-modes.
#' Repeated random initialization (the default 'rnd') appears to work well in most
#' cases.
#'
#' @param data		A data frame of categorical data encoded as non-negative
#'			integers.
#' @param K		Number of clusters to fit.
#' @param algorithm	Algorithm to employ.
#' @param init.method	initialization method to use.
#' @param n.init	Number of random initializations.
#' @param true.column	Column of data frame containing true cluster assignments.
#' @param shuffle	Randomly shuffle observation order before each run.
#' @param verbosity	Level of verbosity.
#' @param seed		Random number seed.
#'
#' @return A list with k-modes clustering solution.
#' \itemize{
#'	\item criterion - The minimum realized objective function value.
#'	\item cluster.criteria - A vector of the minimum realized objective per
#'	cluster.
#'	\item cluster.sizes - The sizes of the clusters in the best solution.
#'	\item partition - The partition of observations into clusters in the best
#'	solution.
#'	\item modes - The estimated modes of each cluster in K x p integer matrix.
#'	\item average.criterion - The average objective function achieved across
#'	multiple initializations.
#'	\item number.initializations - The number of initializations performed.
#'	\item best.ari - If true clusters provided, then the best achieved ARI,
#'	 perhaps not from the best solution, as measured by the objective function.
#'	\item average.ari - If true clusters provided, then the average ARI achieved
#'	across all initializations.
#'	\item rng.seed - The random number seed used for these results, if provided
#'	by the user.
#' }
#'
#' @references
#' \insertRef{Arthur2007}{kmodes}
#' \insertRef{Chaturvedi2001}{kmodes}
#' \insertRef{Cao2009a}{kmodes}
#' \insertRef{devos18}{kmodes}
#' \insertRef{Hartigan1979}{kmodes}
#' \insertRef{Huang1997b}{kmodes}
#' \insertRef{MacQueen1967}{kmodes}
#' \insertRef{Weihs2005}{kmodes}
#'
#' @examples
#' data(sim)
#' sim.h97 <- kmodes(sim, n.init = 10, K = 4, true.column = 1)
#' sim.h97$best.ari
#'
#' ## running Hartigan and Wong algorithm without the true cluster column
#' sim.hw <- kmodes(sim[,-1], algorithm = "hw", n.init = 5, K = 4)
#' if (requireNamespace("mclust", quietly = TRUE)) {
#' 	mclust::adjustedRandIndex(sim[,1], sim.hw$partition)
#' }
#'
#' @export
#' @useDynLib kmodes
kmodes <- function(	data,
			K,
			algorithm = c("h97","hw","lloyd"),
			init.method = c("rnd","h97","h97rnd","hd17","clb09","clb09rnd","av07","av07grd","rndp","rnds"),
			n.init = 1,
			true.column = 0,
			shuffle = TRUE,
			verbosity = 0,
			seed = NULL)
{
	if (is.null(algorithm))
		algorithm <- 'h97'
	if (is.null(init.method))
		init.method <- 'rnd'

	algorithm <- match.arg(algorithm)
	init.method <- match.arg(init.method)

	stopifnot(is.data.frame(data))
	stopifnot(!is.na(as.integer(K)) & as.integer(K)>0)
	stopifnot(is.element(algorithm, c("h97","hw","lloyd")))
	stopifnot(is.element(init.method, c("rnd","h97","h97rnd","hd17","clb09","clb09rnd","av07","av07grd","rndp","rnds")))
	stopifnot(!is.na(as.integer(n.init)) & as.integer(n.init)>0)
	stopifnot(sum(is.na(as.integer(true.column)))==0)
	stopifnot(is.logical(shuffle))
	stopifnot(is.element(as.integer(verbosity), 0:7))

	if (is.null(seed))
		stopifnot(!is.na(as.integer(seed)))

	if (!is.loaded("run_kmodes_r", PACKAGE="kmodes"))
		dyn.load("src/kmodes.so")

	if (true.column == 1) {
		true.cluster <- as.integer(data[, true.column])
		data <- t(as.matrix(data[, (true.column+1):ncol(data)]))
	} else if (true.column == ncol(data)) {
		true.cluster <- as.integer(data[, true.column])
		data <- t(as.matrix(data[, 1:(true.column-1)]))
	} else if (true.column) {
		true.cluster <- as.integer(data[, true.column])
		data <- t(as.matrix(data[, c(1:(true.column-1), (true.column+1):ncol(data))]))
	} else {
		true.cluster <- NULL
		data <- t(as.matrix(data))
	}

	out <- .Call("run_kmodes_r", data, as.integer(K), algorithm, init.method, as.integer(n.init), true.cluster, shuffle, as.integer(verbosity), as.integer(seed))

	return(out)
}# kmodes
