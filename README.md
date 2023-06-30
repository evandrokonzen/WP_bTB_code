# Files to fit model in Konzen *et al.* (2023)

Notation:

* $G$: number of social groups.
* $T$: number of quarters in the whole study, which is on $t = 1, \dots, T$.
* $m$: number of individuals.
* $N$: number of iterations.
* $p$: number of parameters

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

File `finalModel_paper1.R` is used for reading processed data and for choosing priors, initial conditions and details about HMC/RWMH parameter updates.

To fit the model, R package `BIID` (see folder `BIID`) has to be installed. To install this from source
you can use the `devtools` library. (Note that you will need `Rtools` installed on Windows, or the compiler tools on other systems.) For example, if the working directory contains the `BIID` package folder, then running the following should install the package from source.

```
library(devtools)
install("BIID")
```

## Outputs of the fitted model

Outputs (posterior samples for parameters and hidden states) are saved in the directory `resultsDirectory` which is chosen by the user in `finalModel_paper1.R`.

Output files are `[X].Rds`, where `X` is:

* `postPars`: ($N \times p$) matrix with the posterior samples for all model parameters.
* `infTimes`, `infectivityTimes`, `deathTimes`: ($m \times N$) matrices with the infection times, times of infectiousness, and death times, respectively, for each iteration. The value `-10L` denotes the absence of the event.
* `nSusByGroup`, `nExpByGroup`, `nInfByGroup`:  ($G \times T \times N$) arrays showing the number of individuals in each compartment ($S$, $E$, and $I$, respectively) in each social group, at each time point, for each iteration.


