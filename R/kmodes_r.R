##
# @file run_kmodes.R
#
# Purpose: R wrapper for k-modes algorithms available in kmodes C software.
##

#' Run k-modes algorithms.
#'
#' More information here.
#'
#' @param data		data frame of categorical data encoded as non-negative integers
#' @param K		number of clusters (positive integer): 1, 2, ...
#' @param algorithm	algorithm to use (string): h97|hw|hwo
#' @param init.method	initialization method (string): rnd|h97|h97rnd|hd17|clb09|clb09rnd|av07|av07grd|rndp|rnds
#' @param n.init	number of random initializations (positive integer)
#' @param true.column	column in data containing true cluster assignments (integer): 1, 2, ..., ncol(data)
#' @param verbosity	level of verbosity (integer): 0, 1, ..., 7
#' @param seed		random number seed (integer)
#'
#' @example
#' file <- gzfile("data/sim.txt.gz", "rt")
#' d <- read.table(file, header=F)
#' d.h97 <- kmodes(d, n.init = 10, K = 4, true.column = 1)
#' d.h97$best.ari
#' if (!require("mclust", quietly = T)) { install.packages("mclust") }
#' library(mclust)
#' adjustedRandIndex(d[,1], d.h97$partition)
#' d.hw <- kmodes(d, algorithm = "hw", n.init = 10, K = 4, true.column = 1)
#' d.hw$best.ari
#'
#' @export
#' @useDynLib kmodes
kmodes <- function(data, K, algorithm = "h97", init.method = "rnd", n.init = 1, true.column = 0, verbosity = 0, seed = NULL)
{
	stopifnot(is.data.frame(data))
	stopifnot(!is.na(as.integer(K)))
	stopifnot(is.element(algorithm, c("h97","hw","lloyd")))
	stopifnot(is.element(init.method, c("rnd","h97","h97rnd","hd17","clb09","clb09rnd","av07","av07grd","rndp","rnds")))
	stopifnot(!is.na(as.integer(n.init)))
	stopifnot(sum(is.na(as.integer(true.column)))==0)
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

	out <- .Call("run_kmodes_r", data, as.integer(K), algorithm, init.method, as.integer(n.init), true.cluster, as.integer(verbosity), as.integer(seed))

	return(out)
}# kmodes
