## load libraries
library(BIID) # package in folder BIID/
library(tidyverse)
library(reshape2)

## set seed
seed <- 2
set.seed(seed)

## set folder name for results to be stored in
modelFileName <- paste0("outputs_seed", seed, "/")

## set path to data
dataDirectory <- "WPbadgerData/"

#########################################
#####      Set hyperparameters      #####
#########################################

method <- 1  # "HMC"
# method <- 2  # "RWMH"
epsilon <- 0.01
epsilonalphas <- 0.2
epsilonbq <- 0.05
epsilontau <- epsilon # 0.2
epsilonc1 <- 0.05
epsilonsens <- 0.1
lambda_b <- 1
beta_b <- 1
tau_b <- 0.01
a2_b <- 1
b2_b <- 1
c1_b <- 1
xi_mean <- 81
xi_sd <- 60
sd_xi_min <- 8
k <- 1

L <- 30           # only used if method=="HMC"

#########################################
#####      Run MCMC-iFFBS code      #####
#########################################

N <- 25000         # number of MCMC iterations
blockSize <- 1000  # outputs will be saved every 'blockSize' iterations

## Create directory for the outputs (posterior samples and hidden states)
resultsDirectory <- paste0(dataDirectory, modelFileName) # where save Rcpp outputs
if(!dir.exists(resultsDirectory)){dir.create(resultsDirectory)}
methodName <- ifelse(method==1, "HMC", "MH")
path <- paste0(resultsDirectory, methodName, "/")
if(!dir.exists(path)){dir.create(path)}

for(i in seq(1, N, by=blockSize)){
  numFrom <- i
  numTo <- i+blockSize-1
  block <- paste0(path, "Iters_from", sprintf("%.0f", numFrom), 
                  "to", sprintf("%.0f", numTo), "/")
  if(!dir.exists(block)){dir.create(block)}
}

## load data
groupNames <- readRDS(paste0(dataDirectory, "groupNames.rds"))
TestMat <- readRDS(paste0(dataDirectory, "TestMat.rds"))
CaptHist <- readRDS(paste0(dataDirectory, "CaptHist.rds"))
CaptEffort <- readRDS(paste0(dataDirectory, "CaptEffort.rds"))
birthTimes <- readRDS(paste0(dataDirectory, "birthTimes.rds"))
startSamplingPeriod <- readRDS(paste0(dataDirectory, "startSamplingPeriod.rds"))
endSamplingPeriod <- readRDS(paste0(dataDirectory, "endSamplingPeriod.rds"))
capturesAfterMonit <- readRDS(paste0(dataDirectory, "capturesAfterMonit.rds"))
timeVec <- readRDS(paste0(dataDirectory, "timeVec.rds"))
dfTimes <- readRDS(paste0(dataDirectory, "dfTimes.rds"))

## extract summaries from data
G <- length(groupNames)
m <- length(birthTimes)
maxt <- ncol(CaptHist)
testNames <- colnames(TestMat[,-c(1:3)])
numTests <- length(testNames)
numSeasons <- 4

startingQuarter <- 1L  # because timeVec starts on Q1

## Creating Xinit (initial values for the hidden states)
Xinit <- matrix(0L, m, maxt)

# putting NAs before birth date 
for (i in 1:m) {
  if(birthTimes[i]>1){
    Xinit[i, 1:(birthTimes[i]-1)] <- NA
  }
}

# putting infection whenever there's a positive result
for (i in 1:m) {
  df_i <- TestMat[TestMat$idNumber==i, ]
  for (tt in 1:maxt) {
    isInfec <- rep(0L, numTests)
    if(tt %in% df_i$time){
      df_i_tt <- df_i[df_i$time==tt, -c(1:3)]
      
      df_i_alltt <- df_i[, -c(1:3)]
      
      if(any(df_i_tt==1L, na.rm = T)){
        # infection starting at time of first positive test result
        infecTimesEarlier <- 0L
        tStartInfec <- max(birthTimes[i]+1, startSamplingPeriod[i], tt-infecTimesEarlier, 1L)
        Xinit[i, tStartInfec] <- 3L # infection starting at infecTimesEarlier
      }
    }
  }
}

# tauInit <- rgamma(n = 1, shape = 1, rate = 0.1)
tauInit <- runif(n = 1, 4, 16)

# Assuming E becomes I some quarters later and forcing no E->S and I->E
for (i in 1:m) {
  if(any(Xinit[i, ] == 3L, na.rm=T)){
    
    firstInfectedTime <- min(which(Xinit[i, ] == 3L))
    if(firstInfectedTime<maxt){
      Xinit[i, firstInfectedTime:maxt] <- 3L
    }
    
    timefromEtoI <- ceiling(rexp(n = 1, rate = 1/tauInit))
    firstInfectiousTime <- firstInfectedTime + timefromEtoI
    if(firstInfectiousTime<maxt){
      Xinit[i, firstInfectiousTime:maxt] <- 1L
    }
  }
}

# putting 9 after last capture dates
lastCaptureTimes <- sapply(1:m, function(i) {max(which(CaptHist[i,]==1L))})
for(i in 1:m){
  if(lastCaptureTimes[i]<maxt){
    Xinit[ i, (lastCaptureTimes[i]+1):maxt] <- 9L
  }
}

# putting NAs before the start of monitoring period
for (i in 1:m) {
  if(startSamplingPeriod[i]>1L){
    Xinit[i, 1:(startSamplingPeriod[i]-1L)] <- NA
  }
}

# putting NAs after the end of monitoring period
for (i in 1:m) {
  if(endSamplingPeriod[i]<maxt){
    Xinit[i, (endSamplingPeriod[i]+1L):maxt] <- NA
  }
}

for(i in 1:m){
  if((all(is.na(Xinit[i,])))){
    stop("There are individuals with only NAs")
  }
  
  # t0_i <- max(birthTimes[i], 1)
  t0_i <- startSamplingPeriod[i]
  if(any(!(Xinit[i,t0_i:endSamplingPeriod[i]]%in%c(0L,1L,3L,9L)))){
    stop("individuals must have initial values 0,1,3,or 9 during their 
         sampling time period")
  }
}

## Set Brock changepoint
colnames(TestMat)[4] <- "Brock1"
TestMat$Brock2 <- NA # creating a column for Brock2
TestMat <- TestMat[,c(1:4, 9, 5:8)] # reordering columns
testNames <- c("Brock1", "Brock2", testNames[-1])

numTests <- ncol(TestMat) - 3 # updating numTests

# initial Brock changepoint
hp_xi <- c(xi_mean, xi_sd)
xiInit <- min(which(timeVec>=2000))
changePointBrock <- xiInit

## This will be the initial TestMat
for(irow in 1:nrow(TestMat)){
  if(TestMat[irow, "time"]>=changePointBrock){
    if(!is.na(TestMat[irow, "Brock1"])){
      TestMat[irow, "Brock2"] <- TestMat[irow, "Brock1"]
      TestMat[irow, "Brock1"] <- NA
    }
  }
}

# put NAs in test results outside of the monitoring period
count <- 0
for(irow in 1:nrow(TestMat)){
  
  id <- TestMat[irow, "idNumber"]
  time <- TestMat[irow, "time"]
  
  if((time<startSamplingPeriod[id]) || (time>endSamplingPeriod[id])){
    TestMat[irow, c(4:ncol(TestMat))] <- NA
    count <- count + 1
  }
}

# count # number of capture rows where test results will not be used

# checking which animals 
# -- were born before monitoring started in their groups
# -- after 1980

vec <- NULL 
ids <- NULL
for(id in 1:m){
  TestMat_i <- filter(TestMat, idNumber==id)
  time <- TestMat_i[1, "time"] # time of first captured
  g <- TestMat_i[1, "group"] 
  if((birthTimes[id]<startSamplingPeriod[id]) & (birthTimes[id]>1)){
    ids <- c(ids, id)
    vec_id <- c(id, g, birthTimes[id], startSamplingPeriod[id]) 
    vec <- rbind(vec, vec_id)
  }
}
ids
colnames(vec) <- c("id", "firstGroup", "birthTime", "startSamplingPeriod")
rownames(vec) <- NULL
vec <- as.data.frame(vec)

gs <- sort(unique(vec$firstGroup))

starts <- NULL
starts_year <- NULL
for(g_i in seq_along(gs)){
  g <- gs[g_i]
  start_g <- min(which(CaptEffort[g, ]==1))
  start_g_year <- dfTimes$time[dfTimes$idx==start_g]
  starts <- c(starts, start_g)
  starts_year <- c(starts_year, start_g_year)
}

ord <- order(starts)
starts <- starts[ord]
starts_year <- starts_year[ord]
gs <- gs[ord]
nuTimes <- c(1L, as.integer(starts))
numNuTimes <- length(nuTimes)

TestMat_ <- as.matrix(TestMat) # using a matrix to use in Rcpp code
CaptEffort_<- as.matrix(CaptEffort)

thetaNames <- paste0("theta", 1:numTests)
rhoNames <- paste0("rho", 1:numTests)
phiNames <- paste0("phi", 1:numTests)
etaNames <- paste0("eta", 1:numSeasons)
parNames <- c(paste0("alpha",1:G),"lambda", "beta", "q", "tau", 
              "a2", "b2", "c1", 
              paste0("nuE", 0:(numNuTimes-1)),
              paste0("nuI", 0:(numNuTimes-1)),
              "xi",
              thetaNames, rhoNames, phiNames, etaNames)

socGroupSizes <- rep(NA, G)
for(j in 1:G){
  socGroupSizes[j] <- length(unique(filter(TestMat, group==j)$idNumber))
}

K <- median(socGroupSizes)

# hyperparameter values
hp_lambda <- c(1, lambda_b)
hp_beta <- c(1, beta_b)
hp_q <- c(1, 1)
hp_tau <- c(1, tau_b)
hp_a2 <- c(1, a2_b)
hp_b2 <- c(1, b2_b)
hp_c1 <- c(1, c1_b)
hp_nu <- c(1, 1, 1)
# hp_xi initialised above
hp_theta <- c(1, 1)
hp_rho <- c(1, 1)
hp_phi <- c(1, 1)
hp_eta <- c(1, 1)

hp <- list(hp_lambda=hp_lambda, hp_beta=hp_beta, hp_q=hp_q, hp_tau=hp_tau, 
           hp_a2=hp_a2, hp_b2=hp_b2, hp_c1=hp_c1, 
           hp_nu=hp_nu, 
           hp_xi=hp_xi,
           hp_theta=hp_theta, hp_rho=hp_rho, 
           hp_phi=hp_phi, hp_eta=hp_eta)

############################################################################

# Using Rcpp code --------------------------------------

# Choosing initial parameter values from the prior -----------------------

# ## This is equivalent to using 
# initParamValues=Inf
# in which case Rcpp will generates initial values from the prior

lambdaInit <- rgamma(n = 1, shape = 1, rate = 100)
alphaInit <- rgamma(n = 1, shape = 1, rate = 1) # alphastar
betaInit <- rgamma(n = 1, shape = 1, rate = 100)
qInit <- rbeta(n = 1, shape1 = hp_q[1], shape2 = hp_q[2])
# tauInit  # now done above in order to be used in Xinit
a2Init <- rgamma(n = 1, shape = hp_a2[1], rate = 100)
b2Init <- rgamma(n = 1, shape = hp_b2[1], rate = 100)
c1Init <- rgamma(n = 1, shape = hp_c1[1], rate = 100)
nuVecInit <- MCMCpack::rdirichlet(n=numNuTimes, alpha=c(8, 1, 1))
# nuSInit <- nuVecInit[1]
nuEInit <- nuVecInit[,2]
nuIInit <- nuVecInit[,3]
# xiInit was sampled above, and TestMat is constructed given xiInit
thetasInit <- runif(numTests, 0.5, 1)
rhosInit <- runif(numTests, 0.2, 0.8)
phisInit <- runif(numTests, 0.7, 1)
etasInit <- rbeta(n=numSeasons, shape1=hp_eta[1], shape2=hp_eta[2])
initParamValues <- c(alphaInit, lambdaInit, betaInit, qInit, tauInit, 
                     a2Init, b2Init, c1Init, nuEInit, nuIInit, xiInit,
                     thetasInit, rhosInit, phisInit, etasInit)

# initParamValues <- Inf # to generate initial values from the prior

############################################################################

# initParamValues <- Inf # to generate initial values from the prior

## fit model
out_ <- BIID::MCMCiFFBS_(N=N, 
    initParamValues = initParamValues, 
    Xinit=Xinit, 
    TestMat=TestMat_,
    CaptHist=CaptHist, 
    birthTimes=birthTimes,
    startSamplingPeriod=startSamplingPeriod,
    endSamplingPeriod=endSamplingPeriod,
    nuTimes=nuTimes,
    CaptEffort=CaptEffort_,
    capturesAfterMonit=capturesAfterMonit,
    numSeasons=numSeasons, seasonStart=startingQuarter,
    maxt=maxt,
    hp_lambda=hp_lambda, hp_beta=hp_beta, hp_q=hp_q, hp_tau=hp_tau, 
    hp_a2=hp_a2, hp_b2=hp_b2, hp_c1=hp_c1, 
    hp_nu=hp_nu, 
    hp_xi=hp_xi,
    hp_theta=hp_theta, hp_rho=hp_rho,
    hp_phi=hp_phi, hp_eta=hp_eta, k=k, K=K,
    sd_xi_min=sd_xi_min,
    method=method, 
    epsilon=epsilon, 
    epsilonalphas=epsilonalphas, 
    epsilonbq=epsilonbq, 
    epsilontau=epsilontau,
    epsilonc1=epsilonc1, 
    epsilonsens=epsilonsens,
    L=L, 
    path = path, 
    blockSize = blockSize)

colnames(out_) <- parNames
