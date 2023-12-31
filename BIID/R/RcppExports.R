# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

CheckSensSpec_ <- function(numTests, TestField, TestTimes, X) {
    .Call(`_BIID_CheckSensSpec_`, numTests, TestField, TestTimes, X)
}

logS <- function(age, a2, b2, c1) {
    .Call(`_BIID_logS`, age, a2, b2, c1)
}

DlogS_a2 <- function(age, a2, b2) {
    .Call(`_BIID_DlogS_a2`, age, a2, b2)
}

DlogS_b2 <- function(age, a2, b2) {
    .Call(`_BIID_DlogS_b2`, age, a2, b2)
}

DlogS_c1 <- function(age, c1) {
    .Call(`_BIID_DlogS_c1`, age, c1)
}

Dlogpt_a2 <- function(age, a2, b2) {
    .Call(`_BIID_Dlogpt_a2`, age, a2, b2)
}

Dlogpt_b2 <- function(age, a2, b2) {
    .Call(`_BIID_Dlogpt_b2`, age, a2, b2)
}

Dlogpt_c1 <- function(c1) {
    .Call(`_BIID_Dlogpt_c1`, c1)
}

fact <- function(k) {
    .Call(`_BIID_fact`, k)
}

ErlangCDF <- function(x_, k, tau) {
    .Call(`_BIID_ErlangCDF`, x_, k, tau)
}

DerivErlangCDF <- function(x_, k, tau) {
    .Call(`_BIID_DerivErlangCDF`, x_, k, tau)
}

HMC_ <- function(curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, epsilon, L, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) {
    .Call(`_BIID_HMC_`, curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, epsilon, L, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K)
}

HMC_2 <- function(curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, epsilon, epsilonalphas, epsilonbq, epsilontau, epsilonc1, nParsNotGibbs, L, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) {
    .Call(`_BIID_HMC_2`, curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, epsilon, epsilonalphas, epsilonbq, epsilontau, epsilonc1, nParsNotGibbs, L, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K)
}

HMC_thetas_rhos <- function(thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho, epsilon, L) {
    .Call(`_BIID_HMC_thetas_rhos`, thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho, epsilon, L)
}

LocateIndiv <- function(TestMat, birthTimes) {
    .Call(`_BIID_LocateIndiv`, TestMat, birthTimes)
}

#' @title Performs MCMC-iFFBS for CMR-Gompertz-SEID-rho model with group-specific alpha
#' @description Performs MCMC-iFFBS for CMR-Gompertz-SEID-rho model with group-specific alpha
#' @param N Number of MCMC iterations
#' @param Xinit Matrix of dimension (\code{m} x \code{maxt}) containing the initial states
#' @param TestMat Matrix with all capture events. 
#' The columns are (time, id, group, test1, test2, test3, ...). 
#' @param CaptHist Matrix of dimension (\code{m} x \code{maxt}) containing the capture history
#' @param birthTimes Integer vector of length \code{m} representing the birth times
#' @param startSamplingPeriod Integer vector of length \code{m} with the start sampling times
#' @param endSamplingPeriod Integer vector of length \code{m} with the end sampling times
#' @param nuTimes Integer vector saying at which times nu parameters are applied
#' @param CaptEffort (G x maxt) matrix with 1 indicates monitoring and 0 otherwise
#' @param capturesAfterMonit matrix containing the last capture times for individuals captured after the monitoring period
#' @param numSeasons Number of seasons 
#' @param seasonStart Starting season
#' @param maxt Number of time points
#' @param hp_lambda Vector of hyperparameter values for the Gamma prior on mean of alpha (1/lambda)
#' @param hp_beta Vector of hyperparameter values for the Gamma prior on frequency-dependent transmission rate
#' @param hp_q Vector of hyperparameter values for the Gamma prior on q
#' @param hp_tau Vector of hyperparameter values for the Gamma prior on average latent period
#' @param hp_a2 Vector of hyperparameter values for the Gamma prior on Gompertz parameter a2
#' @param hp_b2 Vector of hyperparameter values for the Gamma prior on Gompertz parameter b2
#' @param hp_c1 Vector of hyperparameter values for the Gamma prior on Gompertz parameter c1
#' @param hp_nu Vector of hyperparameter values for Dirichlet prior on the 
#' initial probability of infection of being susceptible, exposed, infectious
#' @param hp_xi Vector of hyperparameter values (mean, std deviation) for the prior on the Brock changepoint
#' @param hp_theta Vector of hyperparameter values for the Beta prior on test 
#' sensitivities
#' @param hp_rho Vector of hyperparameter values for the Beta prior on scaling 
#' factor of test sensitivities in the latent period
#' @param hp_phi Vector of hyperparameter values for the Beta prior on test 
#' specificities
#' @param hp_eta Vector of hyperparameter values for the Beta prior on capture 
#' probabilities
#' @param k Positive integer (shape) parameter of Gamma distribution for the latent period
#' @param K Rescaling parameter for the population size (in order to have beta independent of q)
#' @param sd_xi_min Minimum value for the proposal standard deviation
#' @param method Integer indicating which method is used for updating 
#' infection rates and Gompertz parameters (method=1 for "HMC" and method=2 for "RWMH")
#' @param epsilon Leapfrog stepsize in HMC
#' @param epsilonalphas Leapfrog stepsize in HMC for alphas
#' @param epsilonbq Leapfrog stepsize in HMC for beta and q
#' @param epsilontau Leapfrog stepsize in HMC for tau
#' @param epsilonc1 Leapfrog stepsize in HMC for c1
#' @param epsilonsens Leapfrog stepsize in HMC for thetas and rhos
#' @param L Number of leapfrog steps in HMC
#' @param path Directory's name where results will be saved on
#' @param blockSize number of iterations saved in each folder
#' @param initParamValues Initial parameter values for the two infection rates;
#' rate at latent period; three Gompertz parameters;
#' initial probability of being susceptible, exposed, and infectious;
#' test sensitivities; rho scaling parameters; test specificities; and 
#' seasonal capture probabilities. If set to Inf, then the 
#' initial parameter values will be generated from their priors.
#' @return A matrix containing MCMC draws (\code{N} rows) from
#' the posterior distribution of all parameters.//' 
#' @details \itemize{ \item Columns having test results in \code{TestMat} 
#' should have values
#' \itemize{ 
#' \item \code{0}: negative test 
#' \item \code{1}: positive test 
#' \item \code{NA}: test not observed
#' }
#' \item \code{Capture history} should have values
#' \itemize{
#' \item \code{0}: not captured 
#' \item \code{1}: captured
#' }
#' \item \code{Xinit} should have values
#' \itemize{ 
#' \item \code{NA}: not born yet
#' \item \code{0}: susceptible 
#' \item \code{3}: exposed (infected but not infectious)
#' \item \code{1}: infectious
#' \item \code{9}: dead
#' }
#' }
#' @export
MCMCiFFBS_ <- function(N, initParamValues, Xinit, TestMat, CaptHist, birthTimes, startSamplingPeriod, endSamplingPeriod, nuTimes, CaptEffort, capturesAfterMonit, numSeasons, seasonStart, maxt, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, hp_nu, hp_xi, hp_theta, hp_rho, hp_phi, hp_eta, k, K, sd_xi_min, method, epsilon, epsilonalphas, epsilonbq, epsilontau, epsilonc1, epsilonsens, L, path, blockSize) {
    .Call(`_BIID_MCMCiFFBS_`, N, initParamValues, Xinit, TestMat, CaptHist, birthTimes, startSamplingPeriod, endSamplingPeriod, nuTimes, CaptEffort, capturesAfterMonit, numSeasons, seasonStart, maxt, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, hp_nu, hp_xi, hp_theta, hp_rho, hp_phi, hp_eta, k, K, sd_xi_min, method, epsilon, epsilonalphas, epsilonbq, epsilontau, epsilonc1, epsilonsens, L, path, blockSize)
}

MakeSeasonVec_ <- function(numSeasons, seasonStart, maxt) {
    .Call(`_BIID_MakeSeasonVec_`, numSeasons, seasonStart, maxt)
}

ObsProcess_ <- function(corrector, t0, endTime, id, CaptHist, TestMat_i, TestTimes_i, etas, thetas, rhos, phis, seasonVec) {
    invisible(.Call(`_BIID_ObsProcess_`, corrector, t0, endTime, id, CaptHist, TestMat_i, TestTimes_i, etas, thetas, rhos, phis, seasonVec))
}

RWMH_ <- function(can, curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) {
    .Call(`_BIID_RWMH_`, can, curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K)
}

RWMH_thetas_rhos <- function(thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho, Sigma2) {
    .Call(`_BIID_RWMH_thetas_rhos`, thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho, Sigma2)
}

RWMH_xi <- function(can, cur, hp_xi, TestFieldProposal, TestField, TestTimes, thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod) {
    .Call(`_BIID_RWMH_xi`, can, cur, hp_xi, TestFieldProposal, TestField, TestTimes, thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod)
}

TestMatAsField <- function(TestMat, m) {
    .Call(`_BIID_TestMatAsField`, TestMat, m)
}

TestMatAsFieldProposal <- function(TestFieldProposal, TestField, TestTimes, xi, xiCan, m) {
    invisible(.Call(`_BIID_TestMatAsFieldProposal`, TestFieldProposal, TestField, TestTimes, xi, xiCan, m))
}

TestTimesField <- function(TestMat, m) {
    .Call(`_BIID_TestTimesField`, TestMat, m)
}

TrProbDeath_ <- function(age, a2, b2, c1, logar) {
    .Call(`_BIID_TrProbDeath_`, age, a2, b2, c1, logar)
}

TrProbSurvive_ <- function(age, a2, b2, c1, logar) {
    .Call(`_BIID_TrProbSurvive_`, age, a2, b2, c1, logar)
}

gradThetasRhos <- function(thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho) {
    .Call(`_BIID_gradThetasRhos`, thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho)
}

grad_ <- function(logPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) {
    .Call(`_BIID_grad_`, logPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K)
}

iFFBS_ <- function(alpha_js, b, q, tau, k, K, probDyingMat, LogProbDyingMat, LogProbSurvMat, logProbRest, nuTimes, nuEs, nuIs, thetas, rhos, phis, etas, id, birthTime, startTime, endTime, X, seasonVec, TestMat_i, TestTimes_i, CaptHist, corrector, predProb, filtProb, logTransProbRest, numInfecMat, SocGroup, mPerGroup, idVecAll, logProbStoSgivenSorE, logProbStoEgivenSorE, logProbStoSgivenI, logProbStoEgivenI, logProbStoSgivenD, logProbStoEgivenD, logProbEtoE, logProbEtoI, whichRequireUpdate, sumLogCorrector) {
    invisible(.Call(`_BIID_iFFBS_`, alpha_js, b, q, tau, k, K, probDyingMat, LogProbDyingMat, LogProbSurvMat, logProbRest, nuTimes, nuEs, nuIs, thetas, rhos, phis, etas, id, birthTime, startTime, endTime, X, seasonVec, TestMat_i, TestTimes_i, CaptHist, corrector, predProb, filtProb, logTransProbRest, numInfecMat, SocGroup, mPerGroup, idVecAll, logProbStoSgivenSorE, logProbStoEgivenSorE, logProbStoSgivenI, logProbStoEgivenI, logProbStoSgivenD, logProbStoEgivenD, logProbEtoE, logProbEtoI, whichRequireUpdate, sumLogCorrector))
}

iFFBScalcLogProbRest <- function(i, ttt, logProbRest, X, SocGroup, LogProbDyingMat, LogProbSurvMat, logProbStoSgivenSorE, logProbStoEgivenSorE, logProbStoSgivenI, logProbStoEgivenI, logProbStoSgivenD, logProbStoEgivenD, logProbEtoE, logProbEtoI) {
    invisible(.Call(`_BIID_iFFBScalcLogProbRest`, i, ttt, logProbRest, X, SocGroup, LogProbDyingMat, LogProbSurvMat, logProbStoSgivenSorE, logProbStoEgivenSorE, logProbStoSgivenI, logProbStoEgivenI, logProbStoSgivenD, logProbStoEgivenD, logProbEtoE, logProbEtoI))
}

ivecMinus1 <- function(v) {
    .Call(`_BIID_ivecMinus1`, v)
}

vecSeq <- function(numTests) {
    .Call(`_BIID_vecSeq`, numTests)
}

logPostThetasRhos <- function(thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho) {
    .Call(`_BIID_logPostThetasRhos`, thetas, rhos, X, startSamplingPeriod, endSamplingPeriod, TestField, TestTimes, hp_theta, hp_rho)
}

logPostXi <- function(xiMin, xiMax, xi, hp_xi, TestField_, TestTimes, thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod) {
    .Call(`_BIID_logPostXi`, xiMin, xiMax, xi, hp_xi, TestField_, TestTimes, thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod)
}

logPost_ <- function(logPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) {
    .Call(`_BIID_logPost_`, logPars, G, X, totalNumInfec, SocGroup, totalmPerGroup, birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K)
}

logitD <- function(x) {
    .Call(`_BIID_logitD`, x)
}

logisticD <- function(x) {
    .Call(`_BIID_logisticD`, x)
}

logit <- function(x) {
    .Call(`_BIID_logit`, x)
}

logistic <- function(x) {
    .Call(`_BIID_logistic`, x)
}

sumLogJacobian <- function(xtilde) {
    .Call(`_BIID_sumLogJacobian`, xtilde)
}

multrnorm <- function(mu, Sigma) {
    .Call(`_BIID_multrnorm`, mu, Sigma)
}

logdmultrnorm <- function(x, mu, Sigma) {
    .Call(`_BIID_logdmultrnorm`, x, mu, Sigma)
}

randu <- function(n) {
    .Call(`_BIID_randu`, n)
}

normTransProbRest <- function(logProbs) {
    .Call(`_BIID_normTransProbRest`, logProbs)
}

