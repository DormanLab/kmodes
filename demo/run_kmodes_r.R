#!/usr/bin/Rscript

#system("cd src && R CMD SHLIB -o kmodes_r.so *.c")
source("R/kmodes_r.R")

data(sim)

if (file.exists("/dev/urandom")) {
	seed <- as.integer(system("od -vAn -N4 -td4 < /dev/urandom", intern = T))
} else {
	seed <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5)) * .Machine$integer.max * runif(1)
}

n.init.dm21 <- 50	# initialize DM21 half as many times as other algorithms
n.init.others <- 100

### compare performance of H97 and HW
set.seed(seed)
sim.h97 <- kmodes(sim, n.init = n.init.others, K = 4, true.column = 1, verbosity = 1, seed = seed)
cat("Huang (1997) [H97]\n")
str(sim.h97)

set.seed(seed)
sim.cgc01 <- kmodes(sim, algorithm = "lloyd", n.init = n.init.others, K = 4, true.column = 1, verbosity = 1, seed = seed)
cat("Chaturvedi, Green & Carroll (2001) [CGC01]\n")
str(sim.cgc01)

set.seed(seed)
sim.dm21 <- kmodes(sim, algorithm = "hw", n.init = n.init.dm21, K = 4, true.column = 1, verbosity = 1, seed = seed)
cat("Dorman & Maitra (2021) [DM21]\n")
str(sim.dm21)

cat("\nWith just half as many initializations, DM21 reduces objective on average by",
	sim.h97$average.criterion - sim.dm21$average.criterion,
	"and increases average ARI by",
	sim.dm21$average.ari - sim.h97$average.ari, "over H97.\n\n")

### compare timing of H97, HW, and klaR
if (requireNamespace('klaR', quietly = TRUE) & requireNamespace('rbenchmark', quietly = TRUE)) {
	library(rbenchmark, quietly = TRUE)

	# beware, sometimes this fails for as yet unknown reasons
	print(benchmark(
		"dm21" = kmodes::kmodes(sim[,-1], K=4, alg="hw"),
		"cgc01" = kmodes::kmodes(sim[,-1], K=4, alg="lloyd"),
		"h97" = kmodes::kmodes(sim[,-1], K=4, alg="h97"),
		"klaR" = klaR::kmodes(sim[,-1], modes=4), replications=10))

	benchmark(
		"dm21" = kmodes::kmodes(sim[,-1], K=4, n.init = n.init.dm21, alg="hw"),
		"cgc01" = kmodes::kmodes(sim[,-1], K=4, n.init = n.init.others, alg="lloyd"),
		"h97" = kmodes::kmodes(sim[,-1], K=4, n.init = n.init.others, alg="h97"), replications = 10)
}
