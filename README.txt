
Processed data required to reproduce results:

birthTimes.rds
CaptEffort.rds
CaptHist.rds
capturesAfterMonit.rds
endSamplingPeriod.rds
groupNames.rds
startSamplingPeriod.rds
TestMat.rds
dfTimes.rds
timeVec.rds


File fitting_bTB.R is used for reading processed data and for choosing priors, initial conditions and details about HMC/RWMH parameter updates.

To fit the model, R package BIID (see folder BIID) has to be installed.

Outputs (posterior samples for parameters and hidden states) are saved in the directory "resultsDirectory" which is chosen by the user in finalModel_paper1.R.


