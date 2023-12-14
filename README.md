# Files to fit model in Konzen *et al.* (2023)

## Installation

You will need to install the `BIID` package from source. The package depends on the `Rcpp` and `RcppArmadillo` packages, 
which require the installation of the correct C++ compilers. The guidance below is taken from Sections 2.1.1, 2.1.2 or 
2.1.3 [here](https://teuder.github.io/rcpp4everyone_en/020_install.html).

### Windows

Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html).

(Make sure you tick the option to add Rtools to the PATH whilst installing if requested.)

(Tested on R 4.3.1 on Windows 10.)

### Mac

Install Xcode command line tools. Execute the command `xcode-select --install` in a Terminal.

You might also need to install the gfortran libraries from:

[https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg](https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg)

(Not tested.)

### Linux

Install gcc and related packages (you might also need `gcc-fortran` for some of the dependencies).

In Ubuntu Linux, execute the command `sudo apt-get install r-base-dev` in a Terminal.

(Tested on R 4.3.1 on Ubuntu 22.04.)

### Install package

Once the compilers have been installed, then the version in this repository can be installed from source using the `devtools` package in R. That is, install and/or load the `devtools` package and then set your working directory in R to correspond to the parent directory that contains the `BIID` folder contained in the repository. Then run:

```
library(devtools)
install("BIID")
```

Once installed, the package can be loaded as usual using e.g.

```
library(BIID)
```

### Additional packages

Additional packages required can be installed from within R using:

```
install.packages(c("tidyverse", "reshape2", "MCMCpack"))
```

## Model and data

Notation:

* $G$: number of social groups.
* $T$: number of quarters in the whole study, which is on $t = 1, \dots, T$.
* $m$: number of individuals.
* $N$: number of iterations.
* $p$: number of parameters.
* $J$: number of diagnostic tests.

## Data

Processed data required to reproduce results:

* `birthTimes.rds`: integer vector (of length $m$) with the birth times.
* `CaptEffort.rds`: ($G \times T$) matrix indicating if social group $g$ was being monitored at time $t$ (1 if yes, 0 otherwise).
* `CaptHist.rds`: ($m \times T$) matrix indicating if individual $i$ was captured at time $t$ (1 if yes, 0 otherwise).
* `capturesAfterMonit.rds`: matrix of two columns showing the last capture times (second column) for individuals (first column) captured after the monitoring period. This information is used for the inference of mortality rates.
* `dfTimes.rds`: ($T \times 2$) look-up matrix showing the times in calendar years, e.g.
  
  ```
  > dfTimes
      idx    time
  1     1 1980.00
  2     2 1980.25
  3     3 1980.50
  4     4 1980.75
  5     5 1981.00
  ```
  
* `endSamplingPeriod.rds`: integer vector (of length $m$) showing when the individuals stopped being monitored.
* `groupNames.rds`: character vector (of length $G$) with the social group names.
* `startSamplingPeriod.rds`: integer vector (of length $m$) showing when the individuals started being monitored.
* `TestMat.rds`: matrix where each row represents a capture event. The first three columns are time, individual number, and social group. The remaining columns are diagnostic test results (1 if positive, 0 if negative, `NA` if not conclusive or not taken).
* `timeVec.rds`: vector (of length $T$) showing the dates in calendar years e.g.
  
  ```
  [1] 1980.00 1980.25 1980.50 1980.75 1981.00
  ```

## Fitting the model

File `runmodel.R` is used for reading processed data and for choosing priors, initial conditions and details about HMC/RWMH parameter updates.

## Outputs of the fitted model

Outputs (posterior samples for parameters and hidden states) are saved in the directory `resultsDirectory` which is chosen by the user in `runmodel.R`.

Output files are `[X].Rds`, where `X` is:

* `postPars`: ($N \times p$) matrix with the posterior samples for all model parameters.
* `logLik` and `logPost`: ($N \times 1$) vectors with the log-likelihood and log-posterior, respectively, for all iterations.
* `nSus`, `nExp` and `nInf`: ($T \times N$) matrices showing the total number of individuals in each compartment ($S$, $E$, and $I$, respectively) in the whole population, at each time point, for each iteration.

Other outputs are split into blocks of iterations, with each block stored in a different folder e.g. `Iters_from1to1000`. In the bullet points below $N$ refers to the number of iterations in each block, and then within each folder there are:

* `infTimes`, `infectivityTimes`, `deathTimes`: ($m \times N$) matrices with the infection times, times of infectiousness, and death times, respectively, for each iteration. The value `-10L` denotes the absence of the event.
* `nSusTested`, `nExpTested` and `nInfTested`: ($T \times N \times J$) arrays showing the total number of individuals in each compartment which were tested by each diagnostic test.
* `nSusByGroup`, `nExpByGroup`, `nInfByGroup`:  ($G \times T \times N$) arrays showing the number of individuals in each compartment ($S$, $E$, and $I$, respectively) in each social group, at each time point, for each iteration.
* `nSusTestedByGroup`, `nExpTestedByGroup` and `nInfTestedByGroup`: lists of $G$ elements (one for each social group), where each element is a ($T \times N \times J$) array showing the number of individuals, in each compartment, which were tested by each diagnostic test.
* `AcontribPop`: ($N \times 1$) vector with the relative contribution of the background transmission rates for each iteration in the whole population.
* `AcontribGroup`: ($G \times N$) matrix with the relative contribution of the background transmission rates for each iteration in each social group.
* `AcontribPopTime`: ($T \times N$) matrix with the relative contribution of the background transmission rates for each iteration at each time point in the whole population.
* `AcontribGroupTime`: ($G \times T \times N$) matrix with the relative contribution of the background transmission rates for each iteration at each time point in each social group.
