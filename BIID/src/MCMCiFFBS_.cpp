#include <RcppArmadillo.h>
#include "functions.h"
#include <Rcpp/Benchmark/Timer.h>

//' @title Performs MCMC-iFFBS for CMR-Gompertz-SEID-rho model with group-specific alpha
//' @description Performs MCMC-iFFBS for CMR-Gompertz-SEID-rho model with group-specific alpha
//' @param N Number of MCMC iterations
//' @param Xinit Matrix of dimension (\code{m} x \code{maxt}) containing the initial states
//' @param TestMat Matrix with all capture events. 
//' The columns are (time, id, group, test1, test2, test3, ...). 
//' @param CaptHist Matrix of dimension (\code{m} x \code{maxt}) containing the capture history
//' @param birthTimes Integer vector of length \code{m} representing the birth times
//' @param startSamplingPeriod Integer vector of length \code{m} with the start sampling times
//' @param endSamplingPeriod Integer vector of length \code{m} with the end sampling times
//' @param nuTimes Integer vector saying at which times nu parameters are applied
//' @param CaptEffort (G x maxt) matrix with 1 indicates monitoring and 0 otherwise
//' @param capturesAfterMonit matrix containing the last capture times for individuals captured after the monitoring period
//' @param numSeasons Number of seasons 
//' @param seasonStart Starting season
//' @param maxt Number of time points
//' @param hp_lambda Vector of hyperparameter values for the Gamma prior on mean of alpha (1/lambda)
//' @param hp_beta Vector of hyperparameter values for the Gamma prior on frequency-dependent transmission rate
//' @param hp_q Vector of hyperparameter values for the Gamma prior on q
//' @param hp_tau Vector of hyperparameter values for the Gamma prior on average latent period
//' @param hp_a2 Vector of hyperparameter values for the Gamma prior on Gompertz parameter a2
//' @param hp_b2 Vector of hyperparameter values for the Gamma prior on Gompertz parameter b2
//' @param hp_c1 Vector of hyperparameter values for the Gamma prior on Gompertz parameter c1
//' @param hp_nu Vector of hyperparameter values for Dirichlet prior on the 
//' initial probability of infection of being susceptible, exposed, infectious
//' @param hp_xi Vector of hyperparameter values (mean, std deviation) for the prior on the Brock changepoint
//' @param hp_theta Vector of hyperparameter values for the Beta prior on test 
//' sensitivities
//' @param hp_rho Vector of hyperparameter values for the Beta prior on scaling 
//' factor of test sensitivities in the latent period
//' @param hp_phi Vector of hyperparameter values for the Beta prior on test 
//' specificities
//' @param hp_eta Vector of hyperparameter values for the Beta prior on capture 
//' probabilities
//' @param k Positive integer (shape) parameter of Gamma distribution for the latent period
//' @param K Rescaling parameter for the population size (in order to have beta independent of q)
//' @param sd_xi_min Minimum value for the proposal standard deviation
//' @param method Integer indicating which method is used for updating 
//' infection rates and Gompertz parameters (method=1 for "HMC" and method=2 for "RWMH")
//' @param epsilon Leapfrog stepsize in HMC
//' @param epsilonalphas Leapfrog stepsize in HMC for alphas
//' @param epsilonbq Leapfrog stepsize in HMC for beta and q
//' @param epsilontau Leapfrog stepsize in HMC for tau
//' @param epsilonc1 Leapfrog stepsize in HMC for c1
//' @param epsilonsens Leapfrog stepsize in HMC for thetas and rhos
//' @param L Number of leapfrog steps in HMC
//' @param path Directory's name where results will be saved on
//' @param blockSize number of iterations saved in each folder
//' @param initParamValues Initial parameter values for the two infection rates;
//' rate at latent period; three Gompertz parameters;
//' initial probability of being susceptible, exposed, and infectious;
//' test sensitivities; rho scaling parameters; test specificities; and 
//' seasonal capture probabilities. If set to Inf, then the 
//' initial parameter values will be generated from their priors.
//' @return A matrix containing MCMC draws (\code{N} rows) from
//' the posterior distribution of all parameters.//' 
//' @details \itemize{ \item Columns having test results in \code{TestMat} 
//' should have values
//' \itemize{ 
//' \item \code{0}: negative test 
//' \item \code{1}: positive test 
//' \item \code{NA}: test not observed
//' }
//' \item \code{Capture history} should have values
//' \itemize{
//' \item \code{0}: not captured 
//' \item \code{1}: captured
//' }
//' \item \code{Xinit} should have values
//' \itemize{ 
//' \item \code{NA}: not born yet
//' \item \code{0}: susceptible 
//' \item \code{3}: exposed (infected but not infectious)
//' \item \code{1}: infectious
//' \item \code{9}: dead
//' }
//' }
//' @export
// [[Rcpp::export]]
arma::mat MCMCiFFBS_(int N, 
                     arma::vec initParamValues,
                     const arma::imat& Xinit, 
                     const arma::imat& TestMat, 
                     const arma::imat& CaptHist, 
                     arma::ivec birthTimes, 
                     arma::ivec startSamplingPeriod,
                     arma::ivec endSamplingPeriod,
                     arma::ivec nuTimes,
                     arma::imat CaptEffort,
                     arma::imat capturesAfterMonit,
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
                     arma::vec hp_xi,
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
                     double epsilonsens,
                     int L, 
                     CharacterVector path, 
                     int blockSize) {

  // //
  // Rcpp::Timer timer;
  // int timer_cnt = 0;
  // double prev_time = 0.0;
  // timer.step("start");
  // //

  //
  //
  
  // calling R functions from within C++
  Environment base("package:base");
  Function saveRDS = base["saveRDS"];

  Environment MCMCpack("package:MCMCpack");
  Function ddir = MCMCpack["ddirichlet"];
  Function rdir = MCMCpack["rdirichlet"];
  
  std::string path_ = Rcpp::as<std::string>(path);
  
  if( !((method!=1L) | (method!=2L)) ){
    Rcout << "Please use either method=1 ('HMC') or method=2 ('RWMH')." << std::endl;
  }
  
  unsigned int m = CaptHist.n_rows;
  int numTests = TestMat.n_cols - 3L;
  int G = max(TestMat.col(2));
  
  if((unsigned int)Xinit.n_rows!=m){
    Rcpp::stop("Xinit and CaptHist must include the same number of individuals.");
  }
  if((unsigned int)birthTimes.n_elem!=m){
    Rcpp::stop("birthTimes and CaptHist must include the same number of individuals.");
  }
  if(CaptHist.n_cols!=Xinit.n_cols){
    Rcpp::stop("CaptHist and Xinit must include the same number of time points.");
  }
  
  int numNuTimes = nuTimes.n_elem;
  int nParsNotGibbs = G + 4L + 3L + 2L*numNuTimes + 1L; 
  // G alphas, lambda, beta, q, tau, survival rates, init probs, Brock changepoint
  
  int nParsThetasRhos = 2*numTests; // thetas, rhos

  if(hp_lambda.n_elem!=2L){
    Rcout << "hp_lambda must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_beta.n_elem!=2L){
    Rcout << "hp_beta must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_q.n_elem!=2L){
    Rcout << "hp_q must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_tau.n_elem!=2L){
    Rcout << "hp_tau must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_a2.n_elem!=2L){
    Rcout << "hp_a2 must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_b2.n_elem!=2L){
    Rcout << "hp_b2 must be a vector of 2 hyperparameters" << std::endl;
  }
  if(hp_c1.n_elem!=2L){
    Rcout << "hp_c1 must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_nu.n_elem!=3L){
    Rcout << "hp_nu must be a vector of 3 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_xi.n_elem!=2L){
    Rcout << "hp_xi must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_theta.n_elem!=2L){
    Rcout << "hp_theta must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_rho.n_elem!=2L){
    Rcout << "hp_rho must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_phi.n_elem!=2L){
    Rcout << "hp_phi must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }
  if(hp_eta.n_elem!=2L){
    Rcout << "hp_eta must be a vector of 2 hyperparameters" << std::endl;
    Rcpp::stop(" ");
  }

  NumericVector hp_nu_NumVec = {hp_nu[0], hp_nu[1], hp_nu[2]};
  
  double lambdaInit;
  double alphaStarInit;
  double betaInit;
  double qInit;
  double tauInit;
  double a2Init;
  double b2Init;
  double c1Init;
  arma::vec nuEInit = arma::zeros<arma::vec>(numNuTimes);
  arma::vec nuIInit = arma::zeros<arma::vec>(numNuTimes);
  int xiInit;
  arma::vec thetasInit(numTests);
  arma::vec rhosInit(numTests);
  arma::vec phisInit(numTests);
  arma::vec etasInit(numSeasons);
  
  if(initParamValues.is_finite()){
    
      if((int)(initParamValues.n_elem)!=((nParsNotGibbs-G+1)+3*numTests+numSeasons)){
        Rcout << "If supplied by the user, initParamValues must have " <<
          "correct length. It must be for: " <<  
          "alpha, lambda, beta, q, average latent period, 3 Gompertz parameters, " << 
          "probs of being at E and I compartments at nuTimes, " << 
          "sensitivities, rhos, specificities, and capture probabilities." << std::endl;
        Rcpp::stop(" ");
      }
      Rcout << "Initial parameter values supplied by the user: " << std::endl;
      
      alphaStarInit = initParamValues[0];
      lambdaInit = initParamValues[1];
      betaInit = initParamValues[2];
      qInit = initParamValues[3];
      tauInit = initParamValues[4];
      a2Init = initParamValues[5];
      b2Init = initParamValues[6];
      c1Init = initParamValues[7];
      for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
        nuEInit[i_nu] = initParamValues[7 + 1 + i_nu];
      }
      for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
        nuIInit[i_nu] = initParamValues[7 + 1 + numNuTimes + i_nu];
      }
      xiInit = (int)(initParamValues[7 + 1 + 2*numNuTimes]);
      
      if((xiInit < 1L) || (xiInit > maxt-1L)){
        Rcpp::stop("Initial value for the Brock changepoint is outside the study period.");
      }
      
      for (int iTest=0; iTest<numTests; iTest++) {
        thetasInit[iTest] = initParamValues[nParsNotGibbs-G+1+iTest];
      }
      for (int iTest=0; iTest<numTests; iTest++) {
        rhosInit[iTest] = initParamValues[nParsNotGibbs-G+1+numTests+iTest];
      }
      for (int iTest=0; iTest<numTests; iTest++) {
        phisInit[iTest] = initParamValues[nParsNotGibbs-G+1+2*numTests+iTest];
      }
      for (int s=0; s<numSeasons; s++) {
        etasInit[s] = initParamValues[nParsNotGibbs-G+1+3*numTests+s];
      }
      

      
    
  }else{
      
      Rcout << "Initial parameter values generated from the prior: " << std::endl;
      
      lambdaInit = Rf_rgamma(hp_lambda[0], 1/hp_lambda[1]);
      alphaStarInit = Rf_rgamma(1.0, 1.0);
      betaInit = Rf_rgamma(hp_beta[0], 1/hp_beta[1]);
      qInit = Rf_rbeta(hp_q[0], hp_q[1]);
      tauInit = Rf_rgamma(hp_tau[0], 1/hp_tau[1]);
      a2Init = Rf_rgamma(hp_a2[0], 1/hp_a2[1]);
      b2Init = Rf_rgamma(hp_b2[0], 1/hp_b2[1]);
      c1Init = Rf_rgamma(hp_c1[0], 1/hp_c1[1]);
      
      for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
        NumericVector samp_nuInit = rdir(Named("n")=1L, Named("alpha")=hp_nu_NumVec);
        // nuSInit = samp_nuInit[0];
        nuEInit[i_nu] = samp_nuInit[1];
        nuIInit[i_nu] = samp_nuInit[2];
      }

      xiInit = (int)(Rf_rnorm(hp_xi[0], hp_xi[1]));
      
      int count = 0L;
      while ((xiInit < 1L) || (xiInit > maxt-1L)) {
        count++ ;
        xiInit = (int)(Rf_rnorm(hp_xi[0], hp_xi[1]));
        if(count>100){
          Rcout << "Use a better prior for xi. " <<
            "More than 100 initial values were drawn from the prior, " <<  
            "and none of them was during the study period." << std::endl;
          Rcpp::stop(" ");
        }
      }
      
      
      
      
      for (int iTest=0; iTest<numTests; iTest++) {
        thetasInit[iTest] = Rf_rbeta(hp_theta[0], hp_theta[1]);
      }
      for (int iTest=0; iTest<numTests; iTest++) {
        rhosInit[iTest] = Rf_rbeta(hp_rho[0], hp_rho[1]);
      }
      for (int iTest=0; iTest<numTests; iTest++) {
        phisInit[iTest] = Rf_rbeta(hp_phi[0], hp_phi[1]);
      }
      
      // if(numTests>1L){
      //   phisInit[4] = 1.0;
      // }
      
      for (int s=0; s<numSeasons; s++) {
        etasInit[s] = Rf_rbeta(hp_eta[0], hp_eta[1]);
      }

      

  }

  
  Rcout << "alphaStar = " << alphaStarInit << std::endl;
  Rcout << "lambda = " << lambdaInit << std::endl;
  Rcout << "alpha = " << alphaStarInit*lambdaInit << std::endl;
  Rcout << "beta = " << betaInit << std::endl;
  Rcout << "q = " << qInit << std::endl;
  Rcout << "tau = " << tauInit << std::endl;
  Rcout << "a2 = " << a2Init << std::endl;
  Rcout << "b2 = " << b2Init << std::endl;
  Rcout << "c1 = " << c1Init << std::endl;
  Rcout << "nuEs = " << nuEInit.t() << std::endl;
  Rcout << "nuIs = " << nuIInit.t() << std::endl;
  Rcout << "xi = " << xiInit << std::endl;
  Rcout << "thetas = " << thetasInit.t();
  Rcout << "rhos = " << rhosInit.t();
  Rcout << "phis = " << phisInit.t();
  Rcout << "etas = " << etasInit.t();
  

  // parameter estimates
  int nPars = nParsNotGibbs+3*numTests+numSeasons;
  arma::mat out = arma::zeros<arma::mat>(N, nPars);
  
  // full log-posterior
  arma::vec logPostPerIter = arma::zeros<arma::vec>(N);
  arma::vec logLikPerIter = arma::zeros<arma::vec>(N);
  
  // number of infectives (I), exposed (E), susceptibles (S) and S+E+I at each time point
  arma::imat nSus = arma::zeros<arma::imat>(maxt, N);
  arma::imat nExp = arma::zeros<arma::imat>(maxt, N);
  arma::imat nInf = arma::zeros<arma::imat>(maxt, N);
  

  // number of susceptibles, exposed and infectives which were tested
  arma::icube nSusTested = arma::zeros<arma::icube>(maxt, blockSize, numTests);
  arma::icube nExpTested = arma::zeros<arma::icube>(maxt, blockSize, numTests);
  arma::icube nInfTested = arma::zeros<arma::icube>(maxt, blockSize, numTests);
  
  // now saving these per social group in a field of length G
  arma::icube zerosCube = arma::zeros<arma::icube>(maxt, blockSize, numTests);
  arma::field<arma::icube> nExpTestedPerGroup(G);
  arma::field<arma::icube> nInfTestedPerGroup(G);
  arma::field<arma::icube> nSusTestedPerGroup(G);
  for (int g=0; g<G; g++) {
    nExpTestedPerGroup(g) = zerosCube;
    nInfTestedPerGroup(g) = zerosCube;
    nSusTestedPerGroup(g) = zerosCube;
  }
  
  
  // number of exposed, infectives and total number of individuals per social group 
  // at each time
  arma::icube nSusByGroup = arma::zeros<arma::icube>(G, maxt, blockSize);
  arma::icube nExpByGroup = arma::zeros<arma::icube>(G, maxt, blockSize);
  arma::icube nInfByGroup = arma::zeros<arma::icube>(G, maxt, blockSize);
  
  // infection, infectivity and death times
  arma::imat infTimes = arma::zeros<arma::imat>(m, blockSize);
  arma::imat infectivityTimes = arma::zeros<arma::imat>(m, blockSize);
  arma::imat deathTimes = arma::zeros<arma::imat>(m, blockSize);
  infTimes.fill(-10L);
  infectivityTimes.fill(-10L);
  deathTimes.fill(-10L);
  
  
  // contributions of background rate of infection at times 1,...,T 
  // (i.e. 0,...,maxt-1 indexes in Rcpp) (although at the moment no infection exist in t=1 (t=0 in Rcpp))
  
  arma::mat AcontribPopTime = arma::zeros<arma::mat>(maxt, blockSize);
  arma::vec AcontribPop = arma::zeros<arma::vec>(blockSize);
  
  arma::cube AcontribGroupTime = arma::zeros<arma::cube>(G, maxt, blockSize);
  arma::mat AcontribGroup = arma::zeros<arma::mat>(G, blockSize);
  
  arma::cube AcontribIndivGroupTime = arma::zeros<arma::cube>(m, G, maxt);
  
  // int sliceToSave = 0L;
  // arma::icube Xicube = arma::zeros<arma::icube>(m, maxt, blockSize);
  
  
  // // proportion of times (out of N iterations) that each individual 
  // // became infected
  // arma::vec propInf = arma::zeros<arma::vec>(m);
  
  
  

  // Rcout << "contruct objects: " << std::endl;
  
  
  // infection rates, latent period and Gompertz parameters
  int nParsBlock1 = nParsNotGibbs - 2L*numNuTimes - 1L;  // nuEs, nuIs, and xi are updated separately
  arma::vec pars = arma::zeros<arma::vec>(nParsBlock1);
  
  arma::vec alpha_js = arma::zeros<arma::vec>(G);
  
  for(int g=0; g<G; g++){
    pars[g] = log(alphaStarInit);
    alpha_js[g] = alphaStarInit*lambdaInit;
  }
  pars[G] = log(lambdaInit);
  pars[G+1] = log(betaInit);
  pars[G+2] = logitD(qInit);
  pars[G+3] = log(tauInit);
  pars[G+4] = log(a2Init);
  pars[G+5] = log(b2Init);
  pars[G+6] = log(c1Init);
  
  double lambda = lambdaInit;
  double beta = betaInit;
  
  // qInit = 0.0;
  double q = qInit;
  double ql = logitD(qInit);
  
  double tau = tauInit;
  double a2 = a2Init;
  double b2 = b2Init;
  double c1 = c1Init;
  
  arma::vec nuEs = nuEInit;
  arma::vec nuIs = nuIInit;
  
  int xi = xiInit;
  
  arma::vec thetas = thetasInit;
  arma::vec rhos = rhosInit;
  
  arma::vec phis = phisInit;
  arma::vec etas = etasInit;
  
  arma::ivec seasonVec = MakeSeasonVec_(numSeasons, seasonStart, maxt);

  arma::imat X = Xinit;

  // check last time each individual was captured
  arma::ivec lastCaptureTimes = arma::zeros<arma::ivec>(m);
  for (unsigned int i=0; i<m; i++) {
    arma::ivec capt_hist_i = (CaptHist.row(i)).t();
    arma::uvec whichCapt = arma::find(capt_hist_i == 1L) + 1L;
    if(whichCapt.n_elem>0){
      lastCaptureTimes[i] = max(whichCapt);
    }else{
      lastCaptureTimes[i] = birthTimes[i];
    }
  }


  // Rcout << "lastCaptureTimes: " << std::endl;
  
  // Covariance matrix of posteriors used in RWMH only
  arma::mat Sigma = arma::zeros<arma::mat>(nParsBlock1, nParsBlock1);
  arma::mat Sigma2 = arma::zeros<arma::mat>(2*numTests, 2*numTests);
  Sigma = Sigma.eye(nParsBlock1,nParsBlock1)*0.1;
  arma::vec can = arma::zeros<arma::vec>(nParsBlock1);
  Sigma2 = Sigma2.eye(2*numTests,2*numTests)*0.01;
  
  double sd_xi = 2.0; // starting with a small proposal variance in the first 100 iterations
  
  // double sd_xi = hp_xi[1]; 
  // proposal variance starts with the prior variance and then adapted every
  // 100 iterations
  int count_accept_sd_xi = 0L;
  // if count_accept_sd_xi achieves 3L, then stop adapting sd_xi
  
  
  
  
  
  arma::vec thetas_rhos = arma::zeros<arma::vec>(2*numTests);
  
  // double delta = 1.0;
  
  // Rcout << "Sigma: " << Sigma << std::endl;
  
  // // calculate matrix of ages
  arma::imat ageMat = arma::zeros<arma::imat>(m, maxt);
  ageMat.fill(-10L);
  for(unsigned int i=0; i<m; i++){
    int mint_i = std::max(1L, birthTimes[i]+0L);
    for(unsigned int tt=(mint_i-1L); tt<maxt; tt++){
      ageMat(i, tt) = tt + 1L - birthTimes[i];
    }
  }

  arma::imat SocGroup = LocateIndiv(TestMat, birthTimes); // (m x maxt) matrix

  // Rcout << "SocGroup: " << std::endl;
  
  
  arma::imat CaptHistUsed = CaptHist;
  for(unsigned int ii=0; ii<m; ii++){
    for(unsigned int tt=0; tt<maxt; tt++){
      int g = SocGroup(ii, tt);
      // Rcout << "g = " << g << std::endl;
      
      if(g==0L || g==-1L){
        CaptHistUsed(ii, tt) = 0L;
      }else{
        CaptHistUsed(ii, tt) = CaptHistUsed(ii, tt)*CaptEffort(g-1, tt);  
      }

    }
  }
  
  
  // Rcout << "CaptHistUsed: " << std::endl;
  
  
  /////////////////////////////////////////////////////////////////
  
  // // correcting SocialGroup location for individuals which 
  // // - were captured in a social group where capture effort stopped 
  // // - were captured somewhere else afterwards
  // arma::ivec id_i = TestMat.col(1);
  // 
  // for(unsigned int ir=0; ir<socGroupTermination.n_rows; ir++){
  //   
  //   int id = socGroupTermination(ir,0);
  //   int g = socGroupTermination(ir,1);
  //   
  //   arma::uvec which_i = arma::find(id_i == id);
  //   arma::imat Tests_i = TestMat.rows(which_i);
  //   arma::ivec times_i = Tests_i.col(0);
  //   arma::ivec groups_i = Tests_i.col(2);
  //   
  //   int row = max(arma::find(groups_i == g));  // works because capture events are in order
  //   int captInterrupTime = times_i[row];
  // 
  //   Rcout << "================================= " << std::endl;
  //   Rcout << "id: " << id << std::endl;
  //   Rcout << "g: " << g << std::endl;
  //   Rcout << "times_i: " << times_i.t() << std::endl;
  //   Rcout << "groups_i: " << groups_i.t() << std::endl;
  //   Rcout << "row: " << row << std::endl;
  //   Rcout << "captInterrupTime: " << captInterrupTime << std::endl;
  // 
  //   int gNew = groups_i[row+1];
  //   
  //   int rowNewSocGroup = min(arma::find((groups_i==gNew) && (times_i>captInterrupTime)));
  //   int tmaxCorrection = times_i[rowNewSocGroup];
  // 
  //   Rcout << "tmaxCorrection: " << tmaxCorrection << std::endl;
  //   
  //   
  //   Rcout << "gNew: " << gNew << std::endl;
  //   
  //   for (int tt=(captInterrupTime); tt<tmaxCorrection; tt++) { // captInterrupTime+1 in R
  //     SocGroup(id-1, tt) = gNew;
  //   }
  //   
  // }
  // 
  // CharacterVector SocGroupName = wrap(path_ + "SocGroup.Rds");
  // saveRDS(SocGroup, Named("file", SocGroupName));
  
  /////////////////////////////////////////////////////////////////
  
  
  // numInfecMat is a (G x (maxt-1)) matrix, where the (g,t)-th element
  // is the number of infectious at group g at time t (except the individual 
  // that is being updated).
  // This is for the coming individual to be updated, i.e. the 1st individual here.
  // Inside the function iFFBS, numInfecMat will be updated.
  arma::imat numInfecMat = arma::zeros<arma::imat>(G, maxt-1);
  for(unsigned int tt=0; tt<(maxt-1L); tt++){
      for(unsigned int ii=1; ii<m; ii++){ // without 1st indiv
        if( X(ii,tt)==1L ){
          int g_i_tt = SocGroup(ii,tt);
          if(g_i_tt!=0L){
            numInfecMat(g_i_tt-1L, tt) += 1L;
          }
        }
      }
  }
  
  // mPerGroup is a (G x (maxt)) matrix, where the (g,t)-th element
  // is the number of individuals at group g at time t (except the individual 
  // that is being updated)
  // This is for the coming individual to be updated, i.e. the 1st individual here.
  arma::imat mPerGroup = arma::zeros<arma::imat>(G, maxt);
  for(unsigned int tt=0; tt<maxt; tt++){
    for(unsigned int ii=1; ii<m; ii++){ // without 1st indiv
      if( (X(ii,tt)==0L) || (X(ii,tt)==1L)  || (X(ii,tt)==3L) ){
        int g_i_tt = SocGroup(ii,tt);
        if(g_i_tt!=0L){
          mPerGroup(g_i_tt-1L, tt) += 1L;
        }
      }
    }
  }
  
  arma::field<arma::imat> TestField = TestMatAsField(TestMat, m);
  arma::field<arma::imat> TestFieldProposal = TestField;
  arma::field<arma::ivec> TestTimes = TestTimesField(TestMat, m);
  
  arma::ivec idVecAll = arma::linspace<arma::ivec>(0,m-1,m);
  arma::mat corrector = arma::zeros<arma::mat>(maxt,4);
  arma::mat predProb = arma::zeros<arma::mat>(maxt,4);
  arma::mat filtProb = arma::zeros<arma::mat>(maxt,4);

  arma::mat corrector_theta_rho = arma::zeros<arma::mat>(maxt,4);
  
  // When updating a new individual jj using iFFBS, some transition probabilities of
  // the remaining m-1 individuals have to be updated
  // Below we calculate which individuals have to be updated at the end of the update
  // of individual jj (that is, probs to be used when updating jj+1)
  arma::field<arma::ivec> whichRequireUpdate((m-1)*(maxt-1));
  // the last one is not needed because when updating the 1st individual we 
  // have to recalculate using new infection rates.
  int count = 0L;
  for (int i=0; i<(int)(m-1); i++) {
    
    unsigned int id = i;
    unsigned int idNext = i+1;
    
    for(unsigned int tt=0; tt<maxt-1; tt++){
        std::vector<int> idx;
        for (unsigned int jj=1L; jj<m; jj++) {
          if(((jj!=id)&&(jj!=idNext))&&
             (SocGroup(jj, tt)!=0)&&
             (SocGroup(jj, tt)==SocGroup(id, tt)||
             (SocGroup(jj, tt)==SocGroup(idNext, tt)))){
            idx.push_back(jj);
          }
        }
        arma::ivec idxArma = as<arma::ivec>(wrap(idx));
        whichRequireUpdate(count) = idxArma;
        count += 1L;
    }
    
  }

  
  
  
  // Rcout << "before iter0: " << std::endl;
  
  // //
  // timer.step("before_iterations");
  // timer_cnt++;
  // NumericVector before_iterations(timer);
  // std::cout << "time diff: (before_iterations): " << (before_iterations[timer_cnt] / 1e9) - prev_time << std::endl;
  // prev_time = (before_iterations[timer_cnt] / 1e9);
  // //
  
  int iterSub = 0;
  
  // Start MCMC iterations -------------------------------------------
  for (int iter=0; iter<N; iter++) {
    

    // Rcout << "iter: " << iter << " out of N=" << N << std::endl;
    
    if( (iter>0) & ((iter+1)%1000 == 0)){
      Rcout << "iter: " << iter+1 << " out of N=" << N << std::endl;
      Rcpp::checkUserInterrupt();
    }

    lambda = exp(pars[G]);
    for(int g=0; g<G; g++){
      alpha_js[g] = exp(pars[g])*lambda;
    }
    beta = exp(pars[G+1]);
    q = logisticD(pars[G+2]);
    ql = pars[G+2];
    tau = exp(pars[G+3]);
    a2 = exp(pars[G+4]);
    b2 = exp(pars[G+5]);
    c1 = exp(pars[G+6]);
    
    ql = logitD(q);
    
    double sumLogCorrector = 0.0;
    
    // update death probabilities conditional on age using new Gompertz parameters
    arma::mat probDyingMat = arma::zeros<arma::mat>(m, maxt);
    arma::mat LogProbDyingMat = arma::zeros<arma::mat>(m, maxt);
    arma::mat LogProbSurvMat = arma::zeros<arma::mat>(m, maxt);
    probDyingMat.fill(-10.0);
    for (unsigned int i=0; i<m; i++) {
      for(unsigned int tt=0; tt<maxt; tt++){
        // if(ageMat(i, tt)>-1L){
        if(ageMat(i, tt)>0L){
          if((tt+1)>(unsigned int)lastCaptureTimes[i]){
            double age_i_tt = (double)(ageMat(i,tt));
            double condProbDeath = TrProbDeath_(age_i_tt, a2, b2, c1, false);
            probDyingMat(i, tt) = condProbDeath;
            LogProbDyingMat(i, tt) = log(condProbDeath);
            // LogProbSurvMat(i, tt) = log(1-condProbDeath);
            LogProbSurvMat(i, tt) = TrProbSurvive_(age_i_tt, a2, b2, c1, true);
            
          }else{
            probDyingMat(i, tt) = 0.0;
            LogProbDyingMat(i, tt) = log(0.0);
            LogProbSurvMat(i, tt) = log(1.0);
          }
        }
      }
    }
    
    
    // Rcout << "pars not gibbs in nat scale: " << exp(pars) << std::endl;
    double logProbEtoE = log(1.0 - ErlangCDF(1, k, tau/((double)k)));
    double logProbEtoI = log(ErlangCDF(1, k, tau/((double)k)));
    

    // update probs from tt to tt+1 using new infection rates
    
    
    // arma::mat logProbStoS = arma::zeros<arma::mat>(G, maxt-1);
    // arma::mat logProbStoE = arma::zeros<arma::mat>(G, maxt-1);
    
    
    arma::mat logProbStoSgivenSorE = arma::zeros<arma::mat>(G, maxt-1);
    arma::mat logProbStoEgivenSorE = arma::zeros<arma::mat>(G, maxt-1);
    arma::mat logProbStoSgivenI = arma::zeros<arma::mat>(G, maxt-1);
    arma::mat logProbStoEgivenI = arma::zeros<arma::mat>(G, maxt-1);
    arma::mat logProbStoSgivenD = arma::zeros<arma::mat>(G, maxt-1);
    arma::mat logProbStoEgivenD = arma::zeros<arma::mat>(G, maxt-1);

    for(unsigned int tt=0; tt<(maxt-1L); tt++){
      for(int g=0; g<G; g++){
        
        int mgt;
        double inf_mgt;
        mgt = mPerGroup(g, tt); // without the 1st individual 
        
        if(SocGroup(0,tt)==g+1L){
          
          // if 1st individual is alive and S or E
          inf_mgt = numInfecMat(g, tt)/(pow((double)(mgt+1.0)/K, q));
          logProbStoSgivenSorE(g, tt) = -alpha_js[g] -beta*inf_mgt;
          logProbStoEgivenSorE(g, tt) = Rf_log1mexp(alpha_js[g] + beta*inf_mgt);
          
          // if 1st individual is alive and I
          inf_mgt = (numInfecMat(g, tt)+1.0)/(pow((double)(mgt+1.0)/K, q));
          logProbStoSgivenI(g, tt) = -alpha_js[g] - beta*inf_mgt;
          logProbStoEgivenI(g, tt) = Rf_log1mexp(alpha_js[g] + beta*inf_mgt);
          
          // if 1st individual is dead
          inf_mgt = (numInfecMat(g, tt))/(pow((double)mgt/K, q));
          logProbStoSgivenD(g, tt) = -alpha_js[g] - beta*inf_mgt;
          logProbStoEgivenD(g, tt) = Rf_log1mexp(alpha_js[g] + beta*inf_mgt);
          
        }else{
          
          // if 1st individual is alive and S or E
          inf_mgt = numInfecMat(g, tt)/(pow((double)mgt/K, q));
          double FOI = alpha_js[g] + beta*inf_mgt;
          double log1mexpFOI = Rf_log1mexp(FOI);
          
          logProbStoSgivenSorE(g, tt) = -FOI;
          logProbStoEgivenSorE(g, tt) = log1mexpFOI;
          
          // if 1st individual is alive and I
          logProbStoSgivenI(g, tt) = -FOI;
          logProbStoEgivenI(g, tt) = log1mexpFOI;
          
          // if 1st individual is dead
          logProbStoSgivenD(g, tt) = -FOI;
          logProbStoEgivenD(g, tt) = log1mexpFOI;
          
        }
        
      }
    }
    
      arma::cube logProbRest = arma::zeros<arma::cube>(maxt-1L, 4L, m);
      for (unsigned int jj=1L; jj<m; jj++) {
        for(unsigned int tt=0L; tt<maxt-1L; tt++){
            // update  logProbRest(tt,_,jj) except 1st individual;
            if((X(jj, tt)==0L)||(X(jj, tt)==1L)||(X(jj, tt)==3L)){
              iFFBScalcLogProbRest(jj, tt, logProbRest, X, SocGroup, 
                                   LogProbDyingMat, LogProbSurvMat, 
                                   logProbStoSgivenSorE, logProbStoEgivenSorE, 
                                   logProbStoSgivenI, logProbStoEgivenI, 
                                   logProbStoSgivenD, logProbStoEgivenD, 
                                   logProbEtoE, logProbEtoI);
            }
        }
      }

      arma::mat logTransProbRest = arma::zeros<arma::mat>(maxt-1L, 4L);
      for (unsigned int jj=1L; jj<m; jj++) {
        logTransProbRest += logProbRest.slice(jj);
      }
      
      // Rcout << "sumLogCorrector: " << sumLogCorrector << std::endl;
      
      // //
      // timer.step("before_iFFBS");
      // timer_cnt++;
      // NumericVector before_iFFBS(timer);
      // std::cout << "time diff: (before_iFFBS): " << (before_iFFBS[timer_cnt] / 1e9) - prev_time << std::endl;
      // prev_time = (before_iFFBS[timer_cnt] / 1e9);
      // //
      
      
      // Rcout << "before iFFBS: " << std::endl;
      
      for (unsigned int jj=0; jj<m; jj++) {

      // Rcout << "iter: " << iter << std::endl;
      // Rcout << "jj: " << jj << std::endl;
      
        // updating X(jj, _)
        iFFBS_(alpha_js, beta, q, tau, k, K,
               probDyingMat,
               LogProbDyingMat, 
               LogProbSurvMat,
               logProbRest,
               nuTimes,
               nuEs, 
               nuIs,
               thetas, 
               rhos,
               phis,
               etas, 
               jj+1L, 
               birthTimes[jj],
               startSamplingPeriod[jj],
               endSamplingPeriod[jj],
               X,
               seasonVec,
               TestField(jj),
               TestTimes(jj),
               CaptHist,
               corrector,
               predProb,
               filtProb,
               logTransProbRest,
               numInfecMat, 
               SocGroup,
               mPerGroup,
               idVecAll,
               logProbStoSgivenSorE, logProbStoEgivenSorE, 
               logProbStoSgivenI, logProbStoEgivenI, 
               logProbStoSgivenD, logProbStoEgivenD, 
               logProbEtoE, logProbEtoI, 
               whichRequireUpdate, 
               sumLogCorrector);

      }
      
      // Rcout << "after iFFBS: " << std::endl;
      // Rcout << "iter: " << iter << std::endl;
      // Rcout << "finish iFFBS " << std::endl;
      

      // //
      // timer.step("after_iFFBS");
      // timer_cnt++;
      // NumericVector after_iFFBS(timer);
      // std::cout << "time diff: (after_iFFBS): " << (after_iFFBS[timer_cnt] / 1e9) - prev_time << std::endl;
      // prev_time = (after_iFFBS[timer_cnt] / 1e9);
      // //
      

    arma::ivec lastObsAliveTimes = arma::zeros<arma::ivec>(m);
    for (unsigned int jj=0; jj<m; jj++) {
      arma::uvec which_deadTimes = arma::find(X.row(jj) == 9L);
      if(which_deadTimes.n_elem>0L){
        lastObsAliveTimes[jj] = min(which_deadTimes) + 1L;
      }else{
        // lastObsAliveTimes[jj] = maxt;
        lastObsAliveTimes[jj] = endSamplingPeriod[jj];
      }
    }

    // numInfec here is for all except the 1-st individual 
    // (that is, the vector updated at the iFFBS for the m-th individual).
    // Thus, we need to take into account the 1st individual
    // arma::ivec totalNumInfec = arma::zeros<arma::ivec>(maxt-1L);
    arma::imat totalNumInfec = numInfecMat;
    for(unsigned int tt=0; tt<maxt-1L; tt++){
      if(X(0,tt)==1L){
        int g = SocGroup(0, tt);
        totalNumInfec(g-1, tt) += 1L;
      }
    }
    // similarly for totalmPerGroup:
    arma::imat totalmPerGroup = mPerGroup;
    for(unsigned int tt=0; tt<maxt; tt++){
      if( (X(0,tt)==0L) || (X(0,tt)==1L)  || (X(0,tt)==3L) ){
        int g = SocGroup(0, tt);
        totalmPerGroup(g-1L, tt) += 1L;
      }
    }

      

    // Updating (a, b, tau) and Gompertz parameters using HMC or RWMH

    // Rcout << "start updating pars via HMC " << std::endl;
    
    if(method==1L){
      
      pars = HMC_2(pars, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                  birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit,
                  ageMat, epsilon, epsilonalphas, epsilonbq, epsilontau, epsilonc1, nParsBlock1, L, 
                  hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K);
      
      
    }else if(method==2L){
      
      if((iter>0) & ((iter+1) % 100 == 0)){
        
        int ir0 = floor(iter*0.1);
        // int ir0 = 0L;
        arma::mat histLogFirstPars = arma::zeros<arma::mat>(iter-ir0, nParsBlock1); 

        for (int ir=0; ir<(iter-ir0); ir++) {
          for (int ic=0; ic<nParsBlock1; ic++) {
            if(ic<G){
              histLogFirstPars(ir,ic) = log( out(ir+ir0,ic)/out(ir+ir0,G) );
            }else if(ic==G+2L){
              histLogFirstPars(ir,ic) = logitD( out(ir+ir0,ic) );
            }else{
              histLogFirstPars(ir,ic) = log( out(ir+ir0,ic) );
            }
          }
        }
        arma::mat postVar = arma::cov(histLogFirstPars);
        
        Sigma = (pow(2.38,2)/(postVar.n_cols))*postVar;
        
      }
      
      can = multrnorm(pars, Sigma);
      
      pars = RWMH_(can, pars, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                   birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit,
                   ageMat,
                   hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K);
      
    }
    
    


    // Rcout << "finish updating pars via HMC " << std::endl;
    
    // //
    // timer.step("after_HMC");
    // timer_cnt++;
    // NumericVector after_HMC(timer);
    // std::cout << "time diff: (after_HMC): " << (after_HMC[timer_cnt] / 1e9) - prev_time << std::endl;
    // prev_time = (after_HMC[timer_cnt] / 1e9);
    // //
    //
    // Rcout << "firstPars: " << firstPars << std::endl;
    
    
    ///////////////////////////////////////////////////////////////////////
    // // Updating nuE and nuI using Gibbs Sampling
    // arma::ivec InitXVec = X.col(0);
    // 
    // double InS = sum( InitXVec == 0L );
    // double InE = sum( InitXVec == 3L );
    // double InI = sum( InitXVec == 1L );
    // 
    // NumericVector nuDirichParams = {InS + hp_nu[0], 
    //                                 InE + hp_nu[1],
    //                                 InI + hp_nu[2]};
    // 
    // NumericVector nuSEI = rdir(Named("n")=1L, Named("alpha")=nuDirichParams);
    // // nuS = nuSEI[0];
    // nuE = nuSEI[1];
    // nuI = nuSEI[2];
    
    
    // Updating nuEs and nuIs using Gibbs Sampling
    arma::ivec numS_atnuTimes = arma::zeros<arma::ivec>(numNuTimes);
    arma::ivec numE_atnuTimes = arma::zeros<arma::ivec>(numNuTimes);
    arma::ivec numI_atnuTimes = arma::zeros<arma::ivec>(numNuTimes);
      
    int startTime = 0L; 
    int i_nu = 0L;
    for (unsigned int i=0; i<m; i++) {
      
      startTime = startSamplingPeriod[i];
                                   
      if(birthTimes[i] < startTime){  // born before monitoring started
        
        i_nu = min(arma::find(nuTimes == startTime)); // min is used to convert to int type

        if(X(i,startTime-1L)==0L){
          numS_atnuTimes[i_nu] += 1L;
        }else if(X(i,startTime-1L)==3L){
          numE_atnuTimes[i_nu] += 1L;
        }else if(X(i,startTime-1L)==1L){
          numI_atnuTimes[i_nu] += 1L;
        }
        
      }
      
    }

    NumericMatrix nuDirichParamsMat(numNuTimes, 3);
    NumericMatrix nuSEIMat(numNuTimes, 3);
    for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
      NumericVector nuDirichParams = {numS_atnuTimes[i_nu] + hp_nu[0], 
                                      numE_atnuTimes[i_nu] + hp_nu[1], 
                                      numI_atnuTimes[i_nu] + hp_nu[2]};
      NumericVector nuSEI = rdir(Named("n")=1L, Named("alpha")=nuDirichParams);
      // nuSs[i_nu] = nuSEI[0];
      nuEs[i_nu] = nuSEI[1];
      nuIs[i_nu] = nuSEI[2];
      
      nuDirichParamsMat(i_nu,_) = nuDirichParams;
      nuSEIMat(i_nu,_) = nuSEI;
    }
    

    //////////////////////////////////////////////////////////////////////
    
    
    
    
    // Rcout << "finish updating nuE, nuI " << std::endl;
    
// #    #     #     #    #     #     #    #     #     #    #     #     #    #    

    if(method==1L){
      // HMC step for thetas and rhos
      thetas_rhos = HMC_thetas_rhos(thetas, rhos, X, startSamplingPeriod, 
                                    endSamplingPeriod, TestField,
                                    TestTimes, hp_theta, hp_rho, epsilonsens, L);

    }else if(method==2L){

      if((iter>0) & ((iter+1) % 100 == 0)){
        int ir0 = floor(iter*0.1);
        arma::mat histThetasRhos = arma::zeros<arma::mat>(iter-ir0, nParsThetasRhos);
        for (int ir=0; ir<(iter-ir0); ir++) {
          for (int ic=0; ic<nParsThetasRhos; ic++) {
            double natScaleValue = out(ir+ir0, nParsNotGibbs+ic);
            histThetasRhos(ir,ic) = log( natScaleValue / (1-natScaleValue) );
          }
        }
        arma::mat postVar = arma::cov(histThetasRhos);
        Sigma2 = (pow(2.38,2)/(postVar.n_cols))*postVar;
      }

     // MH step for thetas and rhos
     thetas_rhos = RWMH_thetas_rhos(thetas, rhos, X, startSamplingPeriod, 
                                    endSamplingPeriod, TestField,
                                    TestTimes, hp_theta, hp_rho, Sigma2);

    }

    for (int iTest=0; iTest<numTests; iTest++) {
      thetas[iTest] = thetas_rhos[iTest];
      rhos[iTest] = thetas_rhos[iTest+numTests];
    }
    
  
    // // Fixing test sensitivities
    // thetas[0] = 0.8;
    // thetas[1] = 0.7;
    // thetas[2] = 0.3;
    // rhos[0] = 1.0;
    // rhos[1] = 1.0;
    // rhos[2] = 1.0;
    

    // Updating test specificities using Gibbs Sampling
    arma::imat sensSpecMatrix = CheckSensSpec_(numTests, TestField, 
                                               TestTimes, X);
    
    for (int iTest=0; iTest<numTests; iTest++) {
      phis[iTest] = Rcpp::rbeta(1L,
                                sensSpecMatrix(2,iTest) + hp_phi[0],
                                sensSpecMatrix(3,iTest) + hp_phi[1])[0];

    }
    
    // // Fixing test specificities
    // phis[iTest] = 0.95;
    
    
    
//   #    #     #    #     #    #     #    #     #    #     #    #     #    #  

    // Updating Brock changepoint  -----------
    
    // proposing new xi
    
    
    // // adapting proposal variance (RR approach)
    // if((iter>0) & ((iter+1) % 100 == 0)){
    //   int ir0 = floor(iter*0.1);
    //   // int ir0 = 0L;
    //   arma::vec histXi = arma::zeros<arma::vec>(iter-ir0); 
    //   for (int ir=0; ir<(iter-ir0); ir++) {
    //     histXi[ir] = out(ir+ir0, G+9);
    //   }
    //   double postVar = arma::var(histXi);
    //   sd_xi = sqrt( (pow(2.38,2)/1.0)*postVar );
    // }
    
    // // adapting proposal variance (0.44 acceptance rate target approach)
    if ( (iter>0) & ((iter+1) % 100 == 0) & (count_accept_sd_xi<3L) & (iter<5000)   )  {

      arma::vec out0 = out.col(G+6+1+2*numNuTimes);
      arma::vec vecSub = out0.elem(arma::linspace<arma::uvec>(iter-99,iter,100));
      arma::vec d = diff( vecSub );

      double ccc = 0.0;
      for (unsigned int di=0; di<d.n_elem; di++) {
        if(d[di]!=0){
          ccc += 1.0;
        }
      }

      double acc = ccc/99;

      if (acc < 0.39) {
        sd_xi = 0.9*sd_xi;
        count_accept_sd_xi = 0L;
      }else if (acc > 0.49) {
        sd_xi = 1.1*sd_xi;
        count_accept_sd_xi = 0L;
      }else{
        count_accept_sd_xi += 1L;
      }

      
      if(sd_xi < sd_xi_min){
        sd_xi = sd_xi_min;
      }
      
      Rcout << "iter: " << iter+1L << std::endl;
      Rcout << "sd_xi: " << sd_xi << std::endl;
      
    }
    

    
    int xiCan = (int)(Rf_rnorm(xi, sd_xi));
    
    
    // Rcpp::Rcout << "xi (current)= " << xi << std::endl;
    // Rcpp::Rcout << "xiCan = " << xiCan << std::endl;
    
    // if xiCan==xi, nothing has do be done
    // if xiCan is outside of the studyperiod, then reject it
    if((xiCan!=xi)&&(xiCan >= 1L)&&(xiCan <= maxt-1L)){
      
      // Rcpp::Rcout << "TestMatAsFieldProposal function is used " << std::endl;

      // Function 'TestMatAsFieldProposal' corrects TestFieldProposal given the 
      // current TestField
      // The argument TestFieldProposal in TestMatAsFieldProposal:
      // -- is passed by reference and exists just to avoid creating again 
      // a big object at every iteration. 
      // -- will always be identical to TestField when using as the argument of 
      // function TestMatAsFieldProposal, but inside the function it will be 
      // corrected if xiCan!=xi
      TestMatAsFieldProposal(TestFieldProposal, TestField, TestTimes, xi, xiCan, m);
      


        // Rcpp::Rcout << "After correction: " << std::endl;
        // Rcpp::Rcout << "TestField.row( 1000-1 ) = " << TestField.row( 1000-1 ) << std::endl;
        // Rcpp::Rcout << "TestFieldProposal.row( 1000-1 ) = " << TestFieldProposal.row( 1000-1 ) << std::endl;
        // if(iter==10){
        //   Rcpp::Rcout << "Has TestFieldProposal for id=10 been updated " << std::endl;  
        //   Rcpp::Rcout << "correctly (outside of function) at each iteration?  " << std::endl;  
        //   Rcpp::Rcout << "If so, remove comments " << std::endl;
        //   Rcpp::stop("");
        // }

      

        // depending on the accept-reject step, either TestField or 
        // TestFieldProposal is updated accordingly in the function RWMH_xi:
        
        xi = RWMH_xi(xiCan, xi, hp_xi, TestFieldProposal, TestField, TestTimes, 
                     thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod);

      
    } 
    

      
      
      
      
//   #    #     #    #     #    #     #    #     #    #     #    #     #    #
    
    // // Updating eta using Gibbs Sampling (old version, ignoring irregular trapping)
    // for (int s=0; s<numSeasons; s++) {
    //   
    //   int sumCaptHist_s = 0L;
    //   for(unsigned int tt=0; tt<maxt; tt++){
    //     if(seasonVec[tt]==s+1L){
    //       sumCaptHist_s += sum(CaptHist.col(tt));
    //     }
    //   }
    // 
    //   int sumXAlive_s = 0;
    //   for (unsigned int ir=0; ir<X.n_rows; ir++) {
    //     for (unsigned int ic=0; ic<X.n_cols; ic++) {
    //       if(seasonVec[ic]==s+1L){
    //         if((X(ir,ic)==0L) || (X(ir,ic)==3L) || (X(ir,ic)==1L)){
    //           sumXAlive_s += 1L;
    //         }
    //       }
    //       
    //     }
    //   }
    //   
    //   double sh1_eta = sumCaptHist_s + hp_eta[0];
    //   double sh2_eta = sumXAlive_s - sumCaptHist_s + hp_eta[1];
    //   double eta_s = Rcpp::rbeta(1L, sh1_eta, sh2_eta)[0];
    //   etas[s] = eta_s;
    // }
    
    
    

    
    
    // Updating eta using Gibbs Sampling considering irregular trapping
    
    

    
    
    for (int s=0; s<numSeasons; s++) {
      
      int sumCaptHist_s = 0L;
      for(unsigned int tt=0; tt<maxt; tt++){
        if(seasonVec[tt]==s+1L){
          sumCaptHist_s += sum(CaptHistUsed.col(tt));
        }
      }
      
      int g;
      int sumXAlive_s = 0;
      for (unsigned int ir=0; ir<X.n_rows; ir++) {
        for (unsigned int ic=0; ic<X.n_cols; ic++) {
          if(seasonVec[ic]==s+1L){
            g = SocGroup(ir, ic);
            if(((X(ir,ic)==0L) || (X(ir,ic)==3L) || (X(ir,ic)==1L)) && CaptEffort(g-1, ic)==1L){
              sumXAlive_s += 1L;
            }
          }
          
        }
      }
      
      double sh1_eta = sumCaptHist_s + hp_eta[0];
      double sh2_eta = sumXAlive_s - sumCaptHist_s + hp_eta[1];
      double eta_s = Rcpp::rbeta(1L, sh1_eta, sh2_eta)[0];
      etas[s] = eta_s;
    }
    
    lambda = exp(pars[G]);
    for(int g=0; g<G; g++){
      alpha_js[g] = exp(pars[g])*lambda;
    }
    beta = exp(pars[G+1]);
    q = logisticD(pars[G+2]);
    ql = pars[G+2];
    tau = exp(pars[G+3]);
    a2 = exp(pars[G+4]);
    b2 = exp(pars[G+5]);
    c1 = exp(pars[G+6]);
    
    
    //
    lambda = exp(pars[G]);
    
    for(int g=0; g<G; g++){
      out(iter, g) = exp(pars[g])*lambda;
    }
    out(iter, G) = lambda;
    
    out(iter, G+1) = exp(pars[G+1]);
    out(iter, G+2) = logisticD(pars[G+2]);
    out(iter, G+3) = exp(pars[G+3]);
    out(iter, G+4) = exp(pars[G+4]);
    out(iter, G+5) = exp(pars[G+5]);
    out(iter, G+6) = exp(pars[G+6]);
    
    // out(iter, G+7) = nuE;
    // out(iter, G+8) = nuI;
    
    for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
      out(iter, G+6+1+i_nu) = nuEs[i_nu];
    }
    for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
      out(iter, G+6+1+numNuTimes+i_nu) = nuIs[i_nu];
    }
    out(iter, G+6+1+2*numNuTimes) = xi;
    
    //
    for (int iTest=0; iTest<numTests; iTest++) {
      out(iter, nParsNotGibbs+iTest) = thetas[iTest];
    }
    for (int iTest=0; iTest<numTests; iTest++) {
      out(iter, nParsNotGibbs+numTests+iTest) = rhos[iTest];
    }
    for (int iTest=0; iTest<numTests; iTest++) {
      out(iter, nParsNotGibbs+2*numTests+iTest) = phis[iTest];
    }
    for (int s=0; s<numSeasons; s++) {
      out(iter, nParsNotGibbs+3*numTests+s) = etas[s];
    }
    
    // Rcout << "  " << std::endl;
    // Rcout << "pars in nat scale: " << std::endl;
    // 
    // Rcout << "alpha = " << out(iter,0) << std::endl;
    // Rcout << "beta = " << out(iter,1) << std::endl;
    // Rcout << "tau = " << out(iter,2) << std::endl;
    // Rcout << "a2 = " << out(iter,3) << std::endl;
    // Rcout << "b2 = " << out(iter,4) << std::endl;
    // Rcout << "c1 = " << out(iter,5) << std::endl;
    // Rcout << "thetas = " << thetas.t();
    // Rcout << "phis = " << phis.t();
    // Rcout << "etas = " << etas.t();
    
    // //
    // timer.step("after_Gibbs");
    // timer_cnt++;
    // NumericVector after_Gibbs(timer);
    // std::cout << "time diff: (after_Gibbs): " << (after_Gibbs[timer_cnt] / 1e9) - prev_time << std::endl;
    // prev_time = (after_Gibbs[timer_cnt] / 1e9);
    // //
    
    // ############################################################
    // Rcout << "Finish gibbs samplers " << std::endl;
    
    // Calculating full log-posterior

    
    // The first term includes:
    // * loglik of transition of hidden states;
    // * loglik of the first hidden state (for animals born before study started);
    // * log-prior of (alpha, beta, a2, b2, c2).
    // Note there’s a P(X1) term in iFFBS paper. 
    // These are 1 for born>t0 and S0 otherwise and is calculated inside logPost()
    logPostPerIter[iter] = logPost_(pars, G, X, totalNumInfec, 
                              SocGroup, totalmPerGroup,
                              birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit,
                              ageMat, 
                              hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K);
    
    
    // Calculate logLikPerIter using logPost - logPriors:
    double a_prior = 0.0;
    for(int g=0; g<G; g++){
      double a = exp(pars[g]);
      a_prior += Rf_dgamma(a, 1.0, 1.0, 1L) + log(a);
    }
    lambda = exp(pars[G]);
    beta = exp(pars[G+1]);
    q = logisticD(pars[G+2]);
    ql = pars[G+2];
    tau = exp(pars[G+3]);
    a2 = exp(pars[G+4]);
    b2 = exp(pars[G+5]);
    c1 = exp(pars[G+6]);
    
    
    // Removing Jacobian terms
    logPostPerIter[iter] -= (sum(pars) - pars[G+2]) + (ql - 2*log(1+exp(ql)));
    // rl is in pars, so discarding from the sum(pars)


    double lambda_prior = Rf_dgamma(lambda, hp_lambda[0], 1.0/hp_lambda[1], 1L);
    double b_prior = Rf_dgamma(beta, hp_beta[0], 1.0/hp_beta[1], 1L);
    double q_prior = Rf_dbeta(q, hp_q[0], hp_q[1], 1L);
    double tau_prior = Rf_dgamma(tau, hp_tau[0], 1.0/hp_tau[1], 1L);
    double a2_prior = Rf_dgamma(a2, hp_a2[0], 1.0/hp_a2[1], 1L);
    double b2_prior = Rf_dgamma(b2, hp_b2[0], 1.0/hp_b2[1], 1L);
    double c1_prior = Rf_dgamma(c1, hp_c1[0], 1.0/hp_c1[1], 1L);
    
    double logprior = a_prior + lambda_prior + b_prior + q_prior + tau_prior +
      a2_prior + b2_prior + c1_prior;
    
    logLikPerIter[iter] = logPostPerIter[iter] - logprior;

    // if(((iter+1L) % 10) == 0){
    //   Rcout << "iter = " << iter << std::endl;
    //   
    //   Rcout << "logPostPerIter[iter] begin = " << logPostPerIter[iter]  << std::endl;
    //   Rcout << "logLikPerIter[iter] begin = " << logLikPerIter[iter] << std::endl;
    //   
    // }
    
    // nuEs, nuIs

    NumericVector nuDirichParamsPrior = {hp_nu[0], 
                                         hp_nu[1],
                                         hp_nu[2]};
    
    for (int i_nu=0; i_nu<numNuTimes; i_nu++) {
     NumericVector dens_nu_post_ = ddir(Named("x")=nuSEIMat(i_nu,_), 
                                        Named("alpha")=nuDirichParamsMat(i_nu,_));
      double dens_nu_post = dens_nu_post_[0];
      logPostPerIter[iter] += log(dens_nu_post);
      
      NumericVector dens_nu_prior_ = ddir(Named("x")=nuSEIMat(i_nu,_), 
                                          Named("alpha")=nuDirichParamsPrior);
      double dens_nu_prior = dens_nu_prior_[0];
      logLikPerIter[iter] += log(dens_nu_post) - log(dens_nu_prior);
    }
    
    
    // Adding the observation process part
    logPostPerIter[iter] += sumLogCorrector;
    logLikPerIter[iter] += sumLogCorrector;
    
    
    // Adding remaining log-priors to the log-posterior
    
    // Brock changepoint prior
    double xiLogPrior = Rf_dnorm4(xi, hp_xi[0], hp_xi[1], 1L);
    // thetas, rhos, phis
    double thetasLogPrior = 0.0;
    double rhosLogPrior = 0.0;
    double phisLogPrior = 0.0;
    for (int iTest=0; iTest<numTests; iTest++) {
      thetasLogPrior += Rf_dbeta(thetas[iTest], hp_theta[0], hp_theta[1], 1L);
      rhosLogPrior += Rf_dbeta(rhos[iTest], hp_rho[0], hp_rho[1], 1L);
      phisLogPrior += Rf_dbeta(phis[iTest], hp_phi[0], hp_phi[1], 1L);
    }
    // etas
    double etasLogPrior = 0.0;
    for (int s=0; s<numSeasons; s++) {
      etasLogPrior += Rf_dbeta(etas[s], hp_eta[0], hp_eta[1], 1L);
    }

    
    logPostPerIter[iter] += xiLogPrior + thetasLogPrior + rhosLogPrior + phisLogPrior + 
      etasLogPrior;
      
      
      // //
      // timer.step("after_fullLogPost");
      // timer_cnt++;
      // NumericVector after_fullLogPost(timer);
      // std::cout << "time diff: (after_fullLogPost): " << (after_fullLogPost[timer_cnt] / 1e9) - prev_time << std::endl;
      // prev_time = (after_fullLogPost[timer_cnt] / 1e9);
      // //
      //
      //   //   //   //   //   //   //   //   //   
      
      // Rcout << "Finish logPostPerIter" << std::endl;
      // saving hidden states
      
      // for(unsigned int i=0; i<m; i++){
      //   for(unsigned int tt=0; tt<maxt; tt++){
      //     Xicube(i, tt, iterSub) = X(i, tt);
      //   }
      // }
      
        
      // if(((iter+1L) % 10) == 0){
      //   // Rcout << "iter = " << iter << std::endl;
      // 
      //   Rcout << "a_prior = " << a_prior << std::endl;
      //   Rcout << "lambda_prior = " << lambda_prior << std::endl;
      //   Rcout << "b_prior = " << b_prior << std::endl;
      //   Rcout << "q_prior = " << q_prior << std::endl;
      //   Rcout << "tau_prior = " << tau_prior << std::endl;
      //   Rcout << "a2_prior = " << a2_prior << std::endl;
      //   Rcout << "b2_prior = " << b2_prior << std::endl;
      //   Rcout << "c1_prior = " << c1_prior << std::endl;
      //   
      //   Rcout << "logprior (for first terms) to be discounted = " << logprior << std::endl;
      //   Rcout << "---- " << std::endl;
      //   
      //   Rcout << "InS = " << InS << std::endl;
      //   Rcout << "InE = " << InE << std::endl;
      //   Rcout << "InI = " << InI << std::endl;
      //   
      //   Rcout << "log(dens_nu_prior) = " << log(dens_nu_prior) << std::endl;
      //   
      //   Rcout << "log(dens_nu_post) = " << log(dens_nu_post) << std::endl;
      // 
      //   Rcout << "xi = " << xi << std::endl;
      //   Rcout << "xiLogPrior = " << xiLogPrior << std::endl;
      //   
      //   Rcout << "Rf_dnorm4(xi, hp_xi[0], hp_xi[1], 0L) = " << 
      //     Rf_dnorm4(xi, hp_xi[0], hp_xi[1], 0L) << std::endl;
      //   
      //   Rcout << "logLikPerIter = " << logLikPerIter[iter] << std::endl;
      //   Rcout << "logPostPerIter = " << logPostPerIter[iter] << std::endl;
      //   
      //   Rcout << "------------------------------ " << std::endl;
      // }
        
        
      
      // saving nInf, nSus, nTot, nSusTested, nExpTested, nInfTested
      int g_i_tt;
      for(unsigned int tt=0; tt<maxt; tt++){

          nSus(tt, iter) = sum(X.col(tt) == 0L);
          nExp(tt, iter) = sum(X.col(tt) == 3L);
          nInf(tt, iter) = sum(X.col(tt) == 1L);

          for(unsigned int i=0; i<m; i++){
            g_i_tt = SocGroup(i, tt);
            if( X(i, tt)==0L ){
              nSusByGroup(g_i_tt-1L, tt, iterSub) += 1L;
            }
            if( X(i, tt)==3L ){
              nExpByGroup(g_i_tt-1L, tt, iterSub) += 1L;
            }
            if( X(i, tt)==1L ){
              nInfByGroup(g_i_tt-1L, tt, iterSub) += 1L;
            }
          }
          
          // // nInfByGroup can be obtained from totalNumInfec at all time 
          // // points except at maxt
          // if(tt<(maxt-1)){
          //   for(int g=0; g<G; g++){
          //     nInfByGroup(g, tt, iter) = totalNumInfec(g, tt);
          //   }
          // }else{
          //   // counting at last time point
          //   for(unsigned int i=0; i<m; i++){
          //     if( (X(i, maxt-1)==0L)||(X(i, maxt-1)==3L)||(X(i, maxt-1)==1L) ){
          //       g_i_tt = SocGroup(i, maxt-1);
          //       if( X(i, maxt-1)==1L ){
          //         nInfByGroup(g_i_tt-1L, maxt-1, iter) += 1L;
          //       }
          //     }
          //   }
          // }

      }
      
      // Rcout << "Finish nInfByGroup" << std::endl;
      
      
      
      int tt;
      bool tested;
      for (unsigned int i=0; i<m; i++) {
        arma::imat Tests_i = TestField(i);
        arma::uvec testTimes_i = arma::conv_to<arma::uvec>::from(TestTimes(i)-1L);
        if(testTimes_i.n_elem>0){
          for(unsigned int tt_i=0; tt_i<testTimes_i.n_elem; tt_i++){
            tt = testTimes_i[tt_i];
            for (int iTest=0; iTest<numTests; iTest++) {
                tested = ((Tests_i(tt_i, iTest)==0L)||(Tests_i(tt_i, iTest)==1L));
                if(tested){
                  g_i_tt = SocGroup(i, tt);
                  if(X(i,tt)==1L){
                    nInfTested(tt, iterSub, iTest) += 1L;
                    (nInfTestedPerGroup(g_i_tt-1L))(tt, iterSub, iTest) += 1L;
                  }else if(X(i,tt)==3L){
                    nExpTested(tt, iterSub, iTest) += 1L;
                    (nExpTestedPerGroup(g_i_tt-1L))(tt, iterSub, iTest) += 1L;
                  }else if(X(i,tt)==0L){
                    nSusTested(tt, iterSub, iTest) += 1L;
                    (nSusTestedPerGroup(g_i_tt-1L))(tt, iterSub, iTest) += 1L;
                  }
                }
            }
          }
        }
      }

      // Rcout << "Finish nSusTestedPerGroup" << std::endl;
      
      // arma::cube AcontribIndivGroupTime = arma::zeros<arma::cube>(m, G, maxt);
      AcontribIndivGroupTime.fill(0.0);
      
      
      int g;
      double totalFOI;
      for (unsigned int i=0; i<m; i++) {
        

        if(any(X.row(i)==3L)){
          // propInf[i] += 1L;

          tt = min(arma::find(X.row(i) == 3L));
          infTimes(i, iterSub) = tt + 1L;
          
          if(tt>0){ // at time interval (t-1, t), we do not know infection rates
            g = SocGroup(i, tt) - 1L;
            totalFOI = alpha_js[g] + beta * totalNumInfec(g, tt-1)/(pow((double)(totalmPerGroup(g, tt-1))/K, q));
            AcontribIndivGroupTime(i, g, tt) = alpha_js[g]/totalFOI;
          }
          
        }
        if(any(X.row(i)==1L)){
          infectivityTimes(i, iterSub) = min(arma::find(X.row(i) == 1L)) + 1L;
        }
        if(any(X.row(i)==9L)){
          deathTimes(i, iterSub) = min(arma::find(X.row(i) == 9L)) + 1L;
        }
      }

      
      // Rcout << "Finish deathTimes" << std::endl;
      
      
      
      // arma::vec AcontribPop = arma::zeros<arma::vec>(blockSize);
      double AcontribPop_sum = 0.0;
      double AcontribPop_count = 0.0;
      for (unsigned int i=0; i<m; i++) {
        for (int g=0; g<G; g++) {
          for (unsigned int tt=0; tt<maxt; tt++) {
            if(AcontribIndivGroupTime(i,g,tt)>0.0){
              AcontribPop_sum += AcontribIndivGroupTime(i,g,tt);
              AcontribPop_count += 1.0;
            }
          }
        }
      }
      AcontribPop[iterSub] = AcontribPop_sum/AcontribPop_count;
      
      
      // arma::mat AcontribPopTime = arma::zeros<arma::mat>(maxt, blockSize);
      for (unsigned int tt=0; tt<maxt; tt++) {
        double AcontribPopTime_sum = 0.0;
        double AcontribPopTime_count = 0.0;
        for (unsigned int i=0; i<m; i++) {
          for (int g=0; g<G; g++) {
              if(AcontribIndivGroupTime(i,g,tt)>0.0){
                AcontribPopTime_sum += AcontribIndivGroupTime(i,g,tt);
                AcontribPopTime_count += 1.0;
              }
          }
        }
        AcontribPopTime(tt, iterSub) = AcontribPopTime_sum/AcontribPopTime_count;
      }
      
      
      // arma::mat AcontribGroup = arma::zeros<arma::mat>(G, blockSize);
      for (int g=0; g<G; g++) {
        double AcontribGroup_sum = 0.0;
        double AcontribGroup_count = 0.0;
        for (unsigned int tt=0; tt<maxt; tt++) {
          for (unsigned int i=0; i<m; i++) {
            if(AcontribIndivGroupTime(i,g,tt)>0.0){
              AcontribGroup_sum += AcontribIndivGroupTime(i,g,tt);
              AcontribGroup_count += 1.0;
            }
          }
        }
        AcontribGroup(g, iterSub) = AcontribGroup_sum/AcontribGroup_count;
      }
      
      
      // arma::cube AcontribGroupTime = arma::zeros<arma::cube>(G, maxt, blockSize);
      for (int g=0; g<G; g++) {
        for (unsigned int tt=0; tt<maxt; tt++) {
          double AcontribGroupTime_sum = 0.0;
          double AcontribGroupTime_count = 0.0;
          for (unsigned int i=0; i<m; i++) {
            if(AcontribIndivGroupTime(i,g,tt)>0.0){
              AcontribGroupTime_sum += AcontribIndivGroupTime(i,g,tt);
              AcontribGroupTime_count += 1.0;
            }
          }
          AcontribGroupTime(g, tt, iterSub) = AcontribGroupTime_sum/AcontribGroupTime_count;
        }
      }
      

      
      // Rcout << "Finish after Calc Summaries" << std::endl;
      

      // timer.step("after_Calc_Summaries");
      // timer_cnt++;
      // NumericVector after_Calc_Summaries(timer);
      // std::cout << "time diff: (after_Calc_Summaries): " << (after_Calc_Summaries[timer_cnt] / 1e9) - prev_time << std::endl;
      // prev_time = (after_Calc_Summaries[timer_cnt] / 1e9);
      // std::cout << "------------------- " << std::endl;
      //
      //  
      //   //   //   //   //   //   //   //   //   
      
      
    // ############################################################
    
    // Saving temporary MCMC results in a file called 'fileName'
    if(((iter+1L) % blockSize) == 0){

      Rcout << "File with the first " << iter+1 <<
        " iterations has been saved." << std::endl;
      
      CharacterVector postParsName = wrap(path_ + "postPars.Rds");
      saveRDS(out, Named("file", postParsName));
      
      CharacterVector logPostName = wrap(path_ + "logPost.Rds");
      saveRDS(logPostPerIter, Named("file", logPostName));
      
      CharacterVector logLikName = wrap(path_ + "logLik.Rds");
      saveRDS(logLikPerIter, Named("file", logLikName));
        
      CharacterVector nSusName = wrap(path_ + "nSus.Rds");
      saveRDS(nSus, Named("file", nSusName));
      CharacterVector nExpName = wrap(path_ + "nExp.Rds");
      saveRDS(nExp, Named("file", nExpName));
      CharacterVector nInfName = wrap(path_ + "nInf.Rds");
      saveRDS(nInf, Named("file", nInfName));


      std::string numFrom = std::to_string(iter+1L-blockSize+1L);
      std::string numTo = std::to_string(iter+1L);
      std::string block_ = "Iters_from" + numFrom + "to" + numTo + "/";

      // CharacterVector XicubeName = wrap(path_ + block_ + "Xicube.Rds");
      // saveRDS(Xicube, Named("file", XicubeName));

      CharacterVector nSusTestedName = wrap(path_ + block_ + "nSusTested.Rds");
      saveRDS(nSusTested, Named("file", nSusTestedName));
      CharacterVector nExpTestedName = wrap(path_ + block_ + "nExpTested.Rds");
      saveRDS(nExpTested, Named("file", nExpTestedName));
      CharacterVector nInfTestedName = wrap(path_ + block_ + "nInfTested.Rds");
      saveRDS(nInfTested, Named("file", nInfTestedName));

      CharacterVector nSusTestedPerGroupName = wrap(path_ + block_ + "nSusTestedPerGroup.Rds");
      saveRDS(nSusTestedPerGroup, Named("file", nSusTestedPerGroupName));
      CharacterVector nInfTestedPerGroupName = wrap(path_ + block_ + "nInfTestedPerGroup.Rds");
      saveRDS(nInfTestedPerGroup, Named("file", nInfTestedPerGroupName));
      CharacterVector nExpTestedPerGroupName = wrap(path_ + block_ + "nExpTestedPerGroup.Rds");
      saveRDS(nExpTestedPerGroup, Named("file", nExpTestedPerGroupName));


      CharacterVector nSusByGroupName = wrap(path_ + block_ + "nSusByGroup.Rds");
      saveRDS(nSusByGroup, Named("file", nSusByGroupName));
      CharacterVector nExpByGroupName = wrap(path_ + block_ + "nExpByGroup.Rds");
      saveRDS(nExpByGroup, Named("file", nExpByGroupName));
      CharacterVector nInfByGroupName = wrap(path_ + block_ + "nInfByGroup.Rds");
      saveRDS(nInfByGroup, Named("file", nInfByGroupName));

      CharacterVector infTimesName = wrap(path_ + block_ + "infTimes.Rds");
      saveRDS(infTimes, Named("file", infTimesName));
      CharacterVector infectivityTimesName = wrap(path_ + block_ + "infectivityTimes.Rds");
      saveRDS(infectivityTimes, Named("file", infectivityTimesName));
      CharacterVector deathTimesName = wrap(path_ + block_ + "deathTimes.Rds");
      saveRDS(deathTimes, Named("file", deathTimesName));

      // 
      
      CharacterVector AcontribPopTimeName = wrap(path_ + block_ + "AcontribPopTime.Rds");
      saveRDS(AcontribPopTime, Named("file", AcontribPopTimeName));
      
      CharacterVector AcontribPopName = wrap(path_ + block_ + "AcontribPop.Rds");
      saveRDS(AcontribPop, Named("file", AcontribPopName));
      
      CharacterVector AcontribGroupTimeName = wrap(path_ + block_ + "AcontribGroupTime.Rds");
      saveRDS(AcontribGroupTime, Named("file", AcontribGroupTimeName));
      
      CharacterVector AcontribGroupName = wrap(path_ + block_ + "AcontribGroup.Rds");
      saveRDS(AcontribGroup, Named("file", AcontribGroupName));
      
    }
    
    
    // // Adjusting epsilon for HMC
    // if ( (method==1L) & (iter>0) & ((iter+1) % 100 == 0) & (iter<burnin))  {
    //   
    //   arma::vec out0 = out.col(1);
    //   arma::vec vecSub = out0.elem(arma::linspace<arma::uvec>(iter-99,iter,100));
    //   arma::vec d = diff( vecSub );
    //   
    //   double ccc = 0.0;
    //   for (unsigned int di=0; di<d.n_elem; di++) {
    //     if(d[di]!=0){
    //       // if( abs(d[di]) > 1e-12  ){
    //       ccc += 1.0;
    //     }
    //   }
    //   
    //   double acc = ccc/99;
    //   double y = 1 + 1000*(acc-0.7)*(acc-0.7)*(acc-0.7);
    //   if (y < 0.9) { epsilon = 0.9*epsilon; }
    //   if (y > 1.1) { epsilon = 1.1*epsilon; }
    //   // Rcout << "acc: " << acc << std::endl;
    //   Rcout << "epsilon: " << epsilon << std::endl;
    //   
    // }
    
    
    
    iterSub += 1L;
    
    if(((iter+1L) % blockSize) == 0){
      
      iterSub = 0L;
      
      nSusTested.fill(0L);
      nExpTested.fill(0L);
      nInfTested.fill(0L);
      
      for (int g=0; g<G; g++) {
        nExpTestedPerGroup(g) = zerosCube;
        nInfTestedPerGroup(g) = zerosCube;
        nSusTestedPerGroup(g) = zerosCube;
      }
      
      nSusByGroup.fill(0L);
      nExpByGroup.fill(0L);
      nInfByGroup.fill(0L);
      
      infTimes.fill(-10L);
      infectivityTimes.fill(-10L);
      deathTimes.fill(-10L);
    }

    
  } // end MCMC iterations

  return(out);
}

// List outList(2);
// outList[0] = out;
// outList[1] = Xstored;
// return(outList);
