#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// calculates the corretion/emission probabilities given by the observation process
void ObsProcess_(arma::mat& corrector,
                 int t0,
                 int endTime,
                 int id,
                 const arma::imat& CaptHist, 
                 arma::imat& TestMat_i,
                 arma::ivec& TestTimes_i,
                 arma::vec& etas, 
                 arma::vec& thetas,
                 arma::vec& rhos,
                 arma::vec& phis, 
                 arma::ivec& seasonVec);

// function to calculate a transition probabilities of a single remaining chain given the status of the 
// current individual being updated in the iFFBS.
void iFFBScalcLogProbRest(int i,
                          int ttt,
                          arma::cube& logProbRest,
                          arma::imat& X,
                          arma::imat& SocGroup,
                          arma::mat& LogProbDyingMat, 
                          arma::mat& LogProbSurvMat,
                          arma::mat& logProbStoSgivenSorE, 
                          arma::mat& logProbStoEgivenSorE, 
                          arma::mat& logProbStoSgivenI, 
                          arma::mat& logProbStoEgivenI, 
                          arma::mat& logProbStoSgivenD, 
                          arma::mat& logProbStoEgivenD, 
                          double& logProbEtoE, 
                          double& logProbEtoI) ;

// function to normalise a vector of probabilities which are given in log-scale
arma::rowvec normTransProbRest(arma::rowvec& logProbs);


// function to update hidden states of a single individual using iFFBS
void iFFBS_(arma::vec& alpha_js, 
            double b, double q, double tau, int k, double K, 
            arma::mat& probDyingMat,
            arma::mat& LogProbDyingMat, 
            arma::mat& LogProbSurvMat,
            arma::cube& logProbRest,
            arma::ivec& nuTimes,
            arma::vec& nuEs,
            arma::vec& nuIs,
            arma::vec& thetas, 
            arma::vec& rhos,
            arma::vec& phis,
            arma::vec& etas, 
            int id,
            int birthTime,
            int startTime,
            int endTime,
            arma::imat& X,
            arma::ivec& seasonVec,
            arma::imat& TestMat_i, 
            arma::ivec& TestTimes_i,
            const arma::imat& CaptHist,
            arma::mat& corrector,
            arma::mat& predProb,
            arma::mat& filtProb,
            arma::mat& logTransProbRest,
            arma::imat& numInfecMat, 
            arma::imat& SocGroup,
            arma::imat& mPerGroup,
            arma::ivec& idVecAll, 
            arma::mat& logProbStoSgivenSorE, 
            arma::mat& logProbStoEgivenSorE, 
            arma::mat& logProbStoSgivenI, 
            arma::mat& logProbStoEgivenI, 
            arma::mat& logProbStoSgivenD, 
            arma::mat& logProbStoEgivenD, 
            double& logProbEtoE, 
            double& logProbEtoI,
            arma::field<arma::ivec>& whichRequireUpdate,
            double& sumLogCorrector);

// function to calculate transition probability of dying conditional on age
double TrProbDeath_(double age, double a2, double b2, 
                    double c1, bool logar);

// function to calculate transition probability of surviving conditional on age
double TrProbSurvive_(double age,double a2, double b2, 
                    double c1, bool logar);


// functions to calculate derivatives of logPost with respect to a Gompertz parameter
double Dlogpt_a2(double age, double a2, double b2);

double Dlogpt_b2(double age, double a2, double b2);

double Dlogpt_c1(double c1);


// logarithm of Gompertz survival function at a given age
double logS(double age, double a2, double b2, double c1);

// functions to calculate derivatives of logS with respect to a Gompertz parameter
double DlogS_a2(double age, double a2, double b2);

double DlogS_b2(double age, double a2, double b2);

double DlogS_c1(double age, double c1);

// function to calculate gradients used in HMC
arma::vec grad_(arma::vec& logPars, int G, 
                arma::imat& X, 
                arma::imat& totalNumInfec,
                arma::imat& SocGroup,
                arma::imat& totalmPerGroup,
                arma::ivec& birthTimes,
                arma::ivec& startSamplingPeriod,
                arma::ivec& lastObsAliveTimes, 
                arma::imat& capturesAfterMonit,
                arma::imat& ageMat, 
                arma::vec& hp_lambda,
                arma::vec& hp_beta,
                arma::vec& hp_q,
                arma::vec& hp_tau,
                arma::vec& hp_a2,
                arma::vec& hp_b2,
                arma::vec& hp_c1,
                int k, 
                double K);

// function to calculate log-posteriors used in HMC
double logPost_(arma::vec& logPars, int G, 
                arma::imat& X, 
                arma::imat& totalNumInfec,
                arma::imat& SocGroup,
                arma::imat& totalmPerGroup,
                arma::ivec& birthTimes,
                arma::ivec& startSamplingPeriod,
                arma::ivec& lastObsAliveTimes, 
                arma::imat& capturesAfterMonit,
                arma::imat& ageMat, 
                arma::vec& hp_lambda,
                arma::vec& hp_beta,
                arma::vec& hp_q,
                arma::vec& hp_tau,
                arma::vec& hp_a2,
                arma::vec& hp_b2,
                arma::vec& hp_c1, 
                int k, 
                double K);

// function to perform HMC
arma::vec HMC_(arma::vec curLogPars, int G, 
               arma::imat& X, 
               arma::imat& totalNumInfec,
               arma::imat& SocGroup,
               arma::imat& totalmPerGroup,
               arma::ivec& birthTimes, 
               arma::ivec& startSamplingPeriod,
               arma::ivec& lastObsAliveTimes,
               arma::imat& capturesAfterMonit,
               arma::imat& ageMat,
               double epsilon, double L, 
               arma::vec& hp_lambda,
               arma::vec& hp_beta,
               arma::vec& hp_q,
               arma::vec& hp_tau,
               arma::vec& hp_a2,
               arma::vec& hp_b2,
               arma::vec& hp_c1, 
               int k, double K) ;

arma::vec HMC_2(arma::vec curLogPars, int G, 
                arma::imat& X, 
                arma::imat& totalNumInfec,
                arma::imat& SocGroup,
                arma::imat& totalmPerGroup,
                arma::ivec& birthTimes, 
                arma::ivec& startSamplingPeriod,
                arma::ivec& lastObsAliveTimes,
                arma::imat& capturesAfterMonit,
                arma::imat& ageMat,
                double epsilon, 
                double epsilonalphas, 
                double epsilonbq, 
                double epsilontau,
                double epsilonc1, 
                int nParsNotGibbs, double L, 
                arma::vec& hp_lambda,
                arma::vec& hp_beta,
                arma::vec& hp_q,
                arma::vec& hp_tau,
                arma::vec& hp_a2,
                arma::vec& hp_b2,
                arma::vec& hp_c1,
                int k, double K) ;
  
// function to find number of:
// positive test results for susceptible individuals;
// negative test results for susceptible individuals;
// positive test results for infectious individuals;
// negative test results for infectious individuals;
arma::imat CheckSensSpec_(int numTests, 
                          arma::field<arma::imat>& TestField, 
                          arma::field<arma::ivec>& TestTimes,
                          arma::imat& X);


// function to perform MCMC with iFFBS
arma::mat MCMCiFFBS_(int N, int burnin, 
                     arma::vec initParamValues,
                     const arma::imat& Xinit, 
                     const arma::imat& TestMat, 
                     const arma::imat& CaptHist, 
                     arma::ivec birthTimes, 
                     arma::ivec startSamplingPeriod,
                     arma::ivec endSamplingPeriod,
                     arma::ivec nuTimes,
                     arma::imat CaptEffort,
                     int numSeasons, int seasonStart, 
                     unsigned int maxt, 
                     arma::vec hp_lambda,
                     arma::vec hp_beta,
                     arma::vec hp_q,
                     arma::vec hp_tau,
                     arma::vec hp_a2,
                     arma::vec hp_b2,
                     arma::vec hp_c1,
                     arma::vec hp_nu,
                     double hp_xi,
                     arma::vec hp_theta,
                     arma::vec hp_rho,
                     arma::vec hp_phi,
                     arma::vec hp_eta,
                     int k,
                     double K,
                     double sd_xi_min,
                     int method, 
                     double epsilon, 
                     double epsilonalphas, 
                     double epsilonbq, 
                     double epsilontau,
                     double epsilonc1, 
                     int L, 
                     CharacterVector path, 
                     int blockSize) ;

// function to generate a single draw from a multivariate normal using R function mgcv::rmvn()
arma::vec multrnorm(arma::vec mu, arma::mat Sigma);

// Random Walk MH for infection rates and Gompertz parameters
arma::vec RWMH_(arma::vec can, 
                arma::vec curLogPars, int G, 
                arma::imat& X, 
                arma::imat& totalNumInfec,
                arma::imat& SocGroup,
                arma::imat& totalmPerGroup,
                arma::ivec& birthTimes,
                arma::ivec& startSamplingPeriod,
                arma::ivec& lastObsAliveTimes, 
                arma::imat& capturesAfterMonit,
                arma::imat& ageMat,
                arma::vec& hp_lambda,
                arma::vec& hp_beta,
                arma::vec& hp_q,
                arma::vec& hp_tau,
                arma::vec& hp_a2,
                arma::vec& hp_b2,
                arma::vec& hp_c1,
                int k, double K) ;

// Make vector of length maxt identifying the seasons at times 1,...,maxt
arma::ivec MakeSeasonVec_(int numSeasons, int seasonStart, unsigned int maxt);

// writes the matrix of tests as an arma::field
arma::field<arma::imat> TestMatAsField(const arma::imat& TestMat, int m);

// correct TestFieldProposal for new candidate value for xi
void TestMatAsFieldProposal(arma::field<arma::imat>& TestFieldProposal, 
                            const arma::field<arma::imat>& TestField,
                            const arma::field<arma::ivec>& TestTimes,
                            int xi, int xiCan, int m);
                                               
// writes test times for each individual as an arma::field
arma::field<arma::ivec> TestTimesField(const arma::imat& TestMat, int m);

// creates an (m x maxt) matrix showing to which group each individual belongs to at each time
arma::imat LocateIndiv(const arma::imat& TestMat, arma::ivec& birthTimes);

// RWMH for (thetas, rhos)
arma::vec RWMH_thetas_rhos(arma::vec& thetas, 
                           arma::vec& rhos, 
                           arma::imat& X,
                           arma::ivec& startSamplingPeriod,
                           arma::ivec& endSamplingPeriod,
                           arma::field<arma::imat>& TestField, 
                           arma::field<arma::ivec>& TestTimes,
                           arma::vec& hp_theta,
                           arma::vec& hp_rho,
                           arma::mat Sigma2);

// HMC for (thetas, rhos)
arma::vec HMC_thetas_rhos(arma::vec& thetas, 
                          arma::vec& rhos, 
                          arma::imat& X,
                          arma::ivec& startSamplingPeriod,
                          arma::ivec& endSamplingPeriod,
                          arma::field<arma::imat>& TestField, 
                          arma::field<arma::ivec>& TestTimes,
                          arma::vec& hp_theta,
                          arma::vec& hp_rho,
                          double epsilon,
                          double L);

// logPost of (thetas, rhos)
double logPostThetasRhos(arma::vec& thetas,
                         arma::vec& rhos,
                         arma::imat& X,
                         arma::ivec& startSamplingPeriod,
                         arma::ivec& endSamplingPeriod,
                         arma::field<arma::imat>& TestField, 
                         arma::field<arma::ivec>& TestTimes, 
                         arma::vec& hp_theta,
                         arma::vec& hp_rho);


// gradient of logPost of (thetas, rhos)
arma::vec gradThetasRhos(arma::vec& thetas,
                         arma::vec& rhos,
                         arma::imat& X,
                         arma::ivec& startSamplingPeriod,
                         arma::ivec& endSamplingPeriod,
                         arma::field<arma::imat>& TestField, 
                         arma::field<arma::ivec>& TestTimes, 
                         arma::vec& hp_theta,
                         arma::vec& hp_rho);

// logPost of xi (Brock changepoint) [only using terms which do not cancel out]
double logPostXi(int xiMin,
                 int xiMax,
                 double xi,
                 arma::vec& hp_xi,
                 arma::field<arma::imat>& TestField_, 
                 arma::field<arma::ivec>& TestTimes,
                 arma::vec& thetas,
                 arma::vec& rhos,
                 arma::vec& phis,
                 arma::imat& X,
                 arma::ivec& startSamplingPeriod,
                 arma::ivec& endSamplingPeriod) ;

// RWMH for xi (Brock changepoint)
int RWMH_xi(int can, 
            int cur, 
            arma::vec& hp_xi,
            arma::field<arma::imat>& TestFieldProposal,
            arma::field<arma::imat>& TestField, 
            arma::field<arma::ivec>& TestTimes, 
            arma::vec& thetas,
            arma::vec& rhos,
            arma::vec& phis,
            arma::imat& X,
            arma::ivec& startSamplingPeriod,
            arma::ivec& endSamplingPeriod)  ;
                  
// logit function applied to double
double logitD(double x);

// logistic function applied to double
double logisticD(double x);

// logit function applied to vector
arma::vec logit(arma::vec& x);

// logistic function applied to vector
arma::vec logistic(arma::vec& x);

// sum of logarithm of Jacobian values
double sumLogJacobian(arma::vec& xtilde);

// factorial of k
double fact(int k);

// Erlang CDF (k is the shape and tau is the scale in R::pgamma)
double ErlangCDF(int x_, int k, double tau);

// Derivative of Erlang CDF with respect to tautilde
double DerivErlangCDF(int x_, int k, double tau);


#endif
