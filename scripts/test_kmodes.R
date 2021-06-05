#!/usr/bin/Rscript

source("R/kmodes.R")
load("data/sim.rda")

if (file.exists("/dev/urandom")) {
	seed <- as.integer(system("od -vAn -N4 -td4 < /dev/urandom", intern = T))
} else {
	seed <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5)) * .Machine$integer.max * runif(1)
}

n.init.dm21 <- 50	# initialize DM21 half as many times as other algorithms
n.init.others <- 100

### compare performance of H97 and HW
set.seed(seed)
sim.h97 <- kmodes(sim, n.init = n.init.others, K = 4, true.column = 1, verbosity = 2)
cat("Huang (1997) [H97]\n")
str(sim.h97)

set.seed(seed)
sim.cgc01 <- kmodes(sim, algorithm = "lloyd", n.init = n.init.others, K = 4, true.column = 1, verbosity = 2)
cat("Chaturvedi, Green & Carroll (2001) [CGC01]\n")
str(sim.cgc01)

set.seed(seed)
sim.dm21 <- kmodes(sim, algorithm = "hw", n.init = n.init.dm21, K = 4, true.column = 1, verbosity = 2)
cat("Dorman & Maitra (2021) [DM21]\n")
str(sim.dm21)

cat("\nWith just half as many initializations, DM21 reduces objective on average by",
	sim.h97$average.criterion - sim.dm21$average.criterion,
	"and increases average ARI by",
	sim.dm21$average.ari - sim.h97$average.ari, "over H97.\n\n")
