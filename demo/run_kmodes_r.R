#!/usr/bin/Rscript

#system("cd src && R CMD SHLIB -o kmodes_r.so *.c")
source("R/kmodes_r.R")

file <- gzfile("data/sim.txt.gz", "rt")
d <- read.table(file, header=F)

if (file.exists("/dev/urandom")) {
	seed <- as.integer(system("od -vAn -N4 -td4 < /dev/urandom", intern = T))
} else {
	seed <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5)) * .Machine$integer.max * runif(1)
}

set.seed(seed)
d.h97 <- kmodes(d, n.init = 100, K = 4, true.column = 1, verbosity = 1, seed = seed)
str(d.h97)

set.seed(seed)
d.hw <- kmodes(d, algorithm = "hw", n.init = 100, K = 4, true.column = 1, verbosity = 1, seed = seed)
str(d.hw)
