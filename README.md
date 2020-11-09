k-modes
=======

This software implements algorithms from [Huang (1997)](http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.6.4718), [Chaturvedi _et al._ (2001)](https://doi.org/10.1007/s00357-001-0004-3), and Dorman and Maitra (2020) to minimize the k-modes objective function.
It also implements several initialization methods and K-selection methods.

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Tutorial](#tutorial)
1. [Input Files](#input)
1. [Output](#output)
1. [Troubleshooting](#troubleshooting)
1. [Command-Line Options](#options)
1. [How to Cite](#cite)
1. [Acknowledgements](#acknowledgements)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />

- kmodes requires [cmake](https://cmake.org) (3.5.0 or higher version) and [gcc](https://gcc.gnu.org) (5.4.0 or higher version).
- kmodes requires Rmath, the [R Standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).  Often, the Rmath library (libRmath.a or libRmath.so for Linux or libRmath.dylib for MacOS) will be installed with R, but not always.  Here are some other locations for the library.
	- r-mathlib on [Ubuntu](https://ubuntu.com/) and [Debian](https://www.debian.org/)
	- libRmath on [Fedora](https://ubuntu.com/), [CentOS](https://centos.org/), [Mageia](https://www.mageia.org/en/), and [Mandriva](https://www.openmandriva.org/)
	- Or if all else fails, you can install the Rmath standalone library from the repository [https://github.com/statslabs/rmath](https://github.com/statslabs/rmath)

# Installation <a name = "installation" />

kmodes has been tested under Linux and MacOS.

1. Clone the repository.

    ```sh
    git clone https://github.com/DormanLab/kmodes.git
    ```

2. Configure the project.

   ```sh
   cd kmodes/src
   cmake .
   ```

3. Compile kmodes.  The executable is called ```run_kmodes```.  It will appear in the ```src``` directory you are currently in.

   ```sh
   make
   ```

# Tutorial <a name = "tutorial" />

For this mini-tutorial, we assume you have compiled the executable ```run_kmodes``` and have copied it into the demo directory, where you will carry out this tutorial.
All the files created in this tutorial are in the `demo` directory.  Beware that you will overwrite them if you repeat these commands.

1. Simulate data consisting of 100 observations with 50 coordinates in 4 categories and relative abundances of 0.1, 0.2, 0.25, and 0.45.
Parameters 8.1 and 2.0 determine the amount of 'time' separating the modes and the observations from their modes.
The simulated data are written to file `sim.dat`, **with the true cluster membership of each observation in the first column**, and the simulation settings to `sim.out`.
We set the random seed to 123 and do not actually run k-modes (`-n 0` means 0 initializations).

```
./run_kmodes --simulate 100 50 4 8.1 2.0 --pi 0.10 0.20 0.25 0.45 -f sim.dat -o sim.out -r 123 -n 0
```

2. Estimate the clustering given K=4 clusters using the Hartigan and Wong algorithm (```-w```), outputing the results in ```sim.hw.out```.
We ask for 100 randomly initializations and set the random number seed to 321.
We inform ```run_kmodes``` that the true cluster memberships are given in the first column (```-c 0```), and this information will be used to compute quantities like the adjusted RAND index.

```
./run_kmodes -f sim.dat -o sim.hw.out -k 4 -w -n 100 -i rnd -r 321 -c 0
```

3. Repeat using Huang's (1997) algorithm, which is the default.  Notice, we are using the exact same random initializations, since we set an identical seed.

```
./run_kmodes -f sim.dat -o sim.h97.out -k 4 -n 100 -i rnd -r 321 -c 0
```

4. Compare of the last line of output on my system shows:
```
tail -n 1 sim.h*.out
==> sim.h97.out <==
1843 3.290000 0.093054 1887.360000 8.116644 0.916890 0.014290 0.929269 0.011482 0.102328 0.016519 0.410692 100
==> sim.hw.out <==
1843 1.530000 0.051875 1868.530000 6.822074 0.947685 0.012844 0.959317 0.009732 0.058391 0.013943 0.365033 100

```
These numbers are the minimum value of objective function, the average and standard deviation of the number of iterations, the average and standard deviation of the objective function, the average and standard deviation of the adjusted RAND index, the average and standard deviation of the normalized mutual information, the average and standard deviation of the normalized variation of information, the total time taken, and the number of initializations.
We can see that Huang's algorithm took more time and obtained substantially higher average values of the objective function, but both methods found the true solution handily (see `maximum ARI` in the output file or watch the terminal output as it scrolls past).

5. See if we can do any better by providing the true partition (```-p```).  (Actually, most initializations already find the true clusters.)

```
tail -n 4 sim.out | head -n 1 > sim.part	# extract true partition from simulation output file
./run_kmodes -f sim.dat -o sim.h97.part.out -p sim.part -k 4 -c 0
```

6. Analyze a real dataset.  We run the Hartigan and Wong algorithm for 1000 initializations for each all K in 1 to 8.
Then, we use `run_kmodes` to select K via 9 different K-selection methods.
This last command uses an R script provided in the `scripts` directory.
If that script does not work for you, the effective number of independent coordinates is 15.47, which you can provide directly as the argument to option `-p`.
Though not shown here, you should probably increase the number of initializations as you increase K since optimization is harder for larger K.
```
for k in `seq 1 8`; do ./run_kmodes -f soybean-small.int.data -o soybean-small.K${k}.out soybean-small.K${k}.ini -k $k -w -n 1000 -i rnd -r $RANDOM -c 0; done
./run_kmodes -f soybean-small.int.data --column 0 -k1 soybean-small.K1.out -k2 soybean-small.K2.out -k3 soybean-small.K3.out -k4 soybean-small.K4.out -k5 soybean-small.K5.out -k6 soybean-small.K6.out -k7 soybean-small.K7.out -p `../scripts/eff_p.Rsrc soybean-small.int.data 0.999`
```
The output is not particularly friendly.  Just pay attention to the penultimate line, which in our case (the results will now vary) looks like:
```
soybean-small.int.data     Maxima: J =  7, rJ =  7, kJ =  7, J2 =  3, rJ2 =  1, kJ2 =  3, KL =  3, rKL =  3, kKL =  3
```
where the K selections made for each method is given.  The reported true K for this dataset is 4.

7. To use the Daneel method, you need the timing information for an initialization method without updates during the first iteration.
```
for k in `seq 1 8`; do ./run_kmodes -f soybean-small.int.data -o soybean-small.K${k}.nup.out soybean-small.K${k}.nup.ini -k $k -w -u -n 1000 -i rnd -r $RANDOM -c 0; done
./run_kmodes -f soybean-small.int.data -k1 soybean-small.K1.ini -k2 soybean-small.K2.ini -k3 soybean-small.K3.ini -k4 soybean-small.K4.ini -k5 soybean-small.K5.ini -k6 soybean-small.K6.ini -k7 soybean-small.K7.ini -k8 soybean-small.K8.ini -n1 soybean-small.K1.nup.ini -n2 soybean-small.K2.nup.ini -n3 soybean-small.K3.nup.ini -n4 soybean-small.K4.nup.ini -n5 soybean-small.K5.nup.ini -n6 soybean-small.K6.nup.ini -n7 soybean-small.K7.nup.ini -n8 soybean-small.K8.nup.ini --column 0
```
It is not clear the method works well for these real data, but try the same operations on the simulation data to see how clear the signal can be.


# Preparing input <a name="input" />

The input is a simple text file with one observation per line, no header or comments, and coordinates separated by spaces.
The coordinate values should be non-negative integers.
```run_kmodes``` will issue a warning if they are not contiguous, but you can safely ignore this warning.
If available, the first column can contain the true cluster memberships, which should be provided as integers between 0 and K-1, inclusive and no values unused.
The first 10 lines of the data file created in simulation in the [tutorial](#tutorial) are shown below.
The first column indicates that the first two observations are both in the third cluster.
```
2 2 2 2 0 0 2 1 1 1 1 1 1 1 2 0 2 2 3 1 2 1 0 2 3 3 2 0 2 0 2 2 2 1 1 1 3 3 1 2 1 3 3 2 0 0 2 1 2 1 3
2 1 0 3 3 3 3 1 1 3 2 1 0 1 2 0 3 3 0 0 2 0 0 3 2 2 2 0 3 2 2 0 1 1 2 0 2 3 0 3 1 0 3 2 0 0 1 1 1 1 3
3 1 3 1 1 0 2 0 3 3 3 0 3 1 0 3 0 3 1 2 3 3 2 2 3 0 2 0 1 2 2 1 3 0 1 2 0 1 1 0 3 2 3 2 3 2 3 0 1 1 0
1 1 0 1 3 1 3 3 2 3 2 3 1 2 2 3 0 3 0 0 2 3 0 1 3 2 3 1 2 0 3 2 2 2 3 0 2 3 2 0 3 2 0 3 2 3 2 0 3 3 1
2 3 1 3 2 0 3 1 2 3 1 1 2 0 1 3 1 1 1 0 2 1 0 3 3 2 2 1 2 0 2 0 2 0 2 2 3 1 0 3 1 0 1 0 1 3 0 0 0 1 3
3 0 3 1 1 2 2 0 2 1 3 1 0 2 3 2 2 3 1 2 2 1 2 3 3 0 2 3 1 2 2 2 3 0 2 3 1 1 3 1 1 1 3 1 2 2 3 1 1 1 0
2 1 0 0 2 0 3 1 1 3 1 1 2 3 2 3 0 3 2 1 2 1 0 3 3 0 2 0 2 0 3 3 3 3 0 1 3 3 0 3 0 3 0 0 2 0 2 3 3 1 3
0 3 1 3 1 1 0 3 0 2 1 0 1 2 0 2 0 1 0 2 2 0 3 1 2 3 3 0 2 3 1 3 0 3 2 2 3 0 1 3 2 1 2 1 2 3 1 3 1 2 1
2 1 0 3 1 0 1 2 1 3 3 1 0 2 1 0 1 1 2 0 1 1 0 3 3 3 3 0 2 0 2 3 3 1 2 0 0 3 2 3 1 2 3 1 1 0 2 1 3 3 3
2 1 0 1 0 0 1 1 2 3 1 2 2 0 2 0 0 0 1 1 2 3 1 3 3 2 2 2 1 0 2 1 3 1 3 0 3 1 0 3 1 2 1 0 2 0 3 3 1 1 3
```


# Output Files <a name = "output" />

## Estimating clusters

```run_kmodes``` outputs results to the screen, but you will probably want it to also write the information to a file.
If you provide two arguments to ```-outfile```, free-form information will be output to the first file and the summary of each initialization will be output to the second file.
If you provide  a single argument, both outputs will go to the same output file.

### Informational output
The following information are output in the first output file.  Most entries are key/value pairs.
- `Seed:` The random number generator seed used.
- `Algorithm:` A string identifying the k-modes algorithm used.
- `Initialization:` A string identifying what kind of initialization method was used.
- `Run for N initializations` a line indicating the number of initializations run.  There will be this many lines of output in the second output file.
- `Number of clusters to estimate:` An integer with the user-selected value of K, number of clusters.
- `Initialization output file:` This key is included if the user provided a second filename to the `--outfile` argument.  This file will contain the initialization output.
- If the `--column` argument was provided, then these additional lines are provided:
	- `Number of true clusters:` An integer indicating the true number of clusters.
	- `True cluster assignments:` The next line contains the true cluster assignments of each observation, same order as in the data file.
	- `True cluster sizes:` The next line contains the true sizes of each cluster.
- `Best optimized criterion:` An integer indicating the best value of the objective function obtained.  The next line contains the contribution of each cluster to the overall objective.  The values on the next line sum to the value on this line.
- `Maximum AR:` If the `--column` argument was provided, then this line indicates the maximum achieved Adjusted RAND index.
- `Best cluster sizes:` This line contains the sizes of the clusters in the best solution.  Notice, the clusters may not be ordered in the same way as the true clusters.
- `Best solution originating seeds:` If random initialization was used, this line indicates the indices of the observations that yielded the best solution.
- `Best solution cluster assignments:` The next line will contain the cluster assignments of all observations in the best-scoring solution.
- `Best modes:`  The next K lines specify the modes of the clusters found for the best-scoring solution.
- The last line is a concise summary of the run.  It contains 13 numbers, which are:
	- The minimum value of the objective function obtained.
	- The average number of initializations taken to converge.
	- The standard deviation in the number of initializations taken to converge.
	- The average value of the objective function obtained.
	- The standard deviation in the value of objective function obtained.
	- The average value of the adjusted RAND index (ARI) if `--column` was provided.
	- The standard deviation in the value of ARI if `--column` was provided.
	- The average value of the normalized mutual information (NMI) if `--column` was provided.
	- The standard deviation in the value of NMI if `--column` was provided.
	- The average value of the normalized variation of information (NVI) if `--column` was provided.
	- The standard deviation in the value of NVI if `--column` was provided.
	- The total time (in seconds) taken.
	- The total number of initializations performed.

### Initialization output
The output for each initialization, if directed to a separate file, is suitable for loading into R or a spreadsheet.
Here is a portion of the initialization output for the tutorial:
```
0 2 873 369 487 114 1843 inf 1.000000 -inf 1.000000 0.000000 0.004131
1 2 873 487 369 114 1843 1843 1.000000 1.000000 1.000000 0.000000 0.008764
2 2 487 369 873 114 1843 1843 1.000000 1.000000 1.000000 0.000000 0.013326
3 2 873 114 487 369 1843 1843 1.000000 1.000000 1.000000 0.000000 0.017918
```
The columns displayed are:
- The initialization number.
- The number of iterations to convergence.
- The cost of each cluster (K=4 in this case).
- The total cost.
- The best cost so far, starts at infinity.
- The adjusted RAND index of the current solution, if `--columns` argument provided.
- The best adjusted RAND index so far, if `--columns` argument provided.
- The normalized mutual information of the current solution, if `--columns` argument provided.
- The normalized variation of information of the current solution, if `--columns` argument provided.
- The total compute time used so far (in seconds). 

## Simulation

When simulating data, there are up to two output files, specified via the `--file` and `--outfile` options.

### Simulated data
The simulated data are output in the file given as argument to `--file`.
The format is as described in [Input Files](#input).
In addition, the first column contains the true cluster membership of the observation on that line.

### Simulation information
If requested via command-line argument `--outfile`, information about the simulation may also be output.
The data in this file are key/value pairs.
The meaning of each entry is provided below:
- `Number of observations:` Integer specifying the number of observations simulated.
- `Number of coordinates:` Integer specifying the number of coordinates in each observation.
- `Number of categories:` Integer specifying the number of categories possible at each coordinate.  (Simulation data does not allow the number of categories to vary by coordinate.)
- `Number of true clusters:`  Integer specifying the number of clusters simulated.  In parentheses it reports how many generated observations.
- `CTMC times:` The times used in the continuous time Markov chain simulator.  The first number is the time separating the modes.  The second time is the time separating the observations from the modes.  The more time separating the modes from the ancestor, the easier the clustering.  The more time separating the observations from their modes, the harder the clustering.
- `CTMC probabilities:` The probability a coordinate is altered in the mode relative to the ancestor, and the probability a coordinate is altered in an observation relative to its mode.
- `Mixing proportions:` The probability of each of K clusters provided on the next line.
- `Simulated modes:`  The true modes provided on the next K lines.
- `Mode pairwise distances:` The Hamming distance between each pair of true modes on the next K-1 lines in upper triangle format.
- `Simulated cluster assignments:` The true cluster assignments are on the next line.
- `Simulated cluster sizes:` The size of the simulated clusters are on the next line.
- `Data written to file:` The name of the file where the simulated data were written.


# Troubleshooting <a name = "troubleshooting" />

1. Please note that if you use ```run_kmodes``` to simulate data (```--simulate``` option), then it will output the true cluster identities in the first column.  It is very important that, when doing clustering, you use the ```--column 0``` option to inform ```run_kmodes``` so that it ***does not*** use this column of data for clustering.  Or remove the column before analysis.

2. ```run_kmodes``` is not very flexible about the data it can process.  It assumes that every category is a non-negative integer.

3. If the categories in your dataset are not contiguous numbers, ```run_kmodes``` will warn you.  You can ignore the warning, or use the second argument to option ```--file``` to modify the data file so the warning goes away.

# Command-Line Options <a name = "options" />

Please run `./run_kmodes --help` for detailed information about all available options.  That help is repeated below.
```
RUN_KMODES

NAME
	run_kmodes - cluster observations with categorical coordinates

SYNOPSIS
	run_kmodes [-r SEED -n N -h97|-l|-w -i INI -o OFILE] -k K -f FILE [OPTIONS]

DESCRIPTION
	run_kmodes clusters observations found in IFILE into K clusters.  It
	randomly initialize N times using initialization method INI after
	setting random number seed SEED, and outputs the results in OFILE.

	Data in FILE should be one observation per line, non-negative integers
	for each coordinate, separated by spaces.  There should be no header
	or comments.

	Three algorithms are implemented: Hartigan and Wong's (-w), Huang's
	(-h97), and Lloyd's (-l).

	Several initialization methods are implemented: see the details below.

	You are required to specify the number of clusters K, but run_kmodes
	can help you select K if you run and store the results for multiple K
	and then use options -kK ... multiple times to specify the output
	files.

	You may simulate data (see -s).  Simulated data written to FILE (-f).
	Simulation conditions written to OFILE (-o) if specified.

OPTIONS
   Data input:  These options provide information about the input data.

	-f,--file FILE [FILE2]
		Data are in FILE (with --simulate, simulated data are written
		to FILE). The data will be optionally rewritten to FILE2, after
		possible modification, such as reindexing categories.  Repeat
		filename twice to overwrite the input file.
	-c, --column N
		Set the column N (0-based) containing the true cluster
		assignments (Default: not set).
	-1
		Subtract 1 from the observation categories (Default: no).

  Data output:  These options control where the output goes.

	-o, --outfile OFILE [OFILE2]
		Set the output filename.  If second argument given, split into
		information into first file, run data into second file
		(Default: none).
	-m, --mode MFILE [MFILE2]
		File with true modes.  Write modified modes to MFILE2 if
		specified.

  Run settings: Specify how k-modes will run.

	-k, --K K
		Set the desired number of clusters K.
	-n N
		Set number of initializations [Default: 1].
	-w, --wong
		Run Hartigan and Wong algorithm (Default: no).
		Can combine with -u.
	-qt, --quick-transfer
		Turn off Hartigan & Wong quick-transfer stage.
	-l, --lloyds
		Run Lloyd's algorithm (Default: no).
		Cannot combine with -u.
	--h97
		Run Huang's algorithm (Default: yes).
		Combine with -u to replicate klaR.
	--hartigan
		Use Hartigan updates with Huang's k-modes (Default: no).
		No effect without --h97.
	--shuffle
		Shuffle the observation order on each initialization (Default: no).
	-r, --random SEED
		Set random number seed (Default: 0).
	--continue
		Continue previous run (Default: no).
	-m, --min OPTIMUM
		Record number of initializations and time to this target optimum.
	-q, quiet
		Silence extra output to stderr (Default: yes).
	-h, --help
		This help.

  Initialization: Specify how k-modes will be initialized.

	-i, --initialization
		Random initialization.  Repeat as needed with following arguments.
	   METHOD
		Set initialization method, one of:
		  rnd       K random observations selected as seeds (Default: yes).
		  h97       Huang's initialization method (Default: no).
		  h97rnd    Huang's initialization method randomized (Default: no).
		  hd17      Huang's initialization method interpreted by de Vos (Default: no).
		  clb09     Cao et al.'s initialization method (Default: no).
		  clb09rnd  Cao et al.'s initialization method randomized (Default: no).
		  av07      k-modes++ (Default: no).
		  av07grd   greedy k-modes++ (Default: no).
		  rndp      given partition via --column, select one
		            observation from each partition (Default: no).
		  rnds      given seed observations, randomly select K;
		            if K seeds provided, this is deterministic (Default: no).
	   INT1 ... INTK
		Set the (0-based) indices of the seeds.
	   IFILE
		Provide file with possible seeds.
		If more than K seeds in IFILE, then method is 'rnds'.
		If K seeds in IFILE, then initialize with these seeds.
	-p, --partition PFILE
		Partition file for deterministic initialization.  Modes of given
		partition will be initial modes.  Partition file contains space-
		separated, 0-based indices of cluster assignments for each
		observation in the order they appear in the input FILE.  This
		option overrides -i.
	-t, --time SECONDS
		Number of seconds to run initializations (Default: 0.00).
	-u, --update
		Use mode updates during first iteration.
		Combine with -i rnd to replicate klaR (Default: no).

  K selection: Perform K-selection based on previous runs.

	-kK FILE1 FILE2 ...
		The names of output files from various runs/methods for K=K,
		limited by POSIX ARG_MAX (excuse unconventional option).
		For DM method, must also use -nK arguments for same values of K.
	-nK FILE1 FILE2 ...
		The names of output files from various runs/methods with NO updates
		during the first iteration for K=K, limited by POSIX_ARG_MAX.
		Must be used with -kK for same values of K.
	-p FLOAT
		Effective number of independent coordinates.

  Simulation: Specify how to simulate data.  Simulation K inferred from -pi.

	-pi, --pi
		Use one of following formats:
	   PROPORTION1 ... PROPORTIONK
		The cluster proportions in simulation.
	   dir FLOAT1 ... FLOATK
		The alpha for Dirichlet prior on pi.  (Type "dir" literally.)
	-s N P C T1 T2
		Simulation size (observations N by coordinates P by categories C)
		and times:
		  T1 is the time separating the "modes".
		  T2 is the time separating the observations from their "modes".

Note to self: Use OFILE2 (-o) and MFILE2 (-m) to convert from fqmorph output
format to run_kmodes format, where categories use contiguous indices 0, 1, ....
Note to self: hidden options --run

RUN_KMODES
```


# How to Cite <a name = "citing" />

- This work is under review.  Please see [arxiv](url to be added).

# Acknowledgements <a name = "acknowledgements" />

- This software makes use of the [R8lib library](https://people.sc.fsu.edu/~jburkardt/c_src/r8lib/r8lib.html), which was released under [GNU LGPL](https://www.gnu.org/licenses/lgpl-3.0.en.html).
- MacQueen's version of the k-modes algorithm was published in [Huang (1997)](http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.6.4718).
- Lloyd's version of the k-modes algorithm was published in [Chaturvedi _et al._ (2001)](https://doi.org/10.1007/s00357-001-0004-3).
- The Hartigan and Wong algorithm for k-means was published in [Hartigan and Wong (1979)](https://www.jstor.org/stable/2346830).
- The h97 initialization method was published in [Huang (1997)](http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.6.4718).
The randomized version h97rnd was suggested by us.
Another randomized version hd17 was suggested by the author of [Python k-modes](https://pypi.org/project/kmodes/)
- The clb09 initialization method was published in [Cao _et al._  (2009)](http://dx.doi.org/10.1016/j.eswa.2009.01.060).  The randomized version was suggested by us.
- The initialization method k-means++ was published in [Arthur and Vassilvitskii (2007)](http://dl.acm.org/citation.cfm?id=1283383.1283494).
A greedy version of k-means++ was proposed by the same authors.
The initialization method k-modes++ and its greedy version were inspired by k-means++, and they were proposed by us.
- The adjusted RAND index was published in [Hubert and Arabie (1985)](http://dx.doi.org/10.1007/BF01908075).
- Normalized mutual information was published in [Kvalseth (1987)](https://ieeexplore.ieee.org/document/4309069).
- Normalized variation of information was published in [Vinh _et al._ (2010)](https://dl.acm.org/doi/10.5555/1756006.1953024).
- MacQueen's and Lloyd's algorithms for k-modes are available in [R package klaR](https://cran.r-project.org/web/packages/klaR/index.html).
- MacQueen's algorithm for k-modes is available in [Python](https://pypi.org/project/kmodes/).


# Contact <a name = "contact" />

If you have any problems with this k-modes software, please contact:

Karin Dorman (kdorman@iastate.edu)
