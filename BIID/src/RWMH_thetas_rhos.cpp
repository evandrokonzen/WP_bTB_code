#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::vec RWMH_thetas_rhos(arma::vec& thetas, 
                           arma::vec& rhos, 
                           arma::imat& X,
                           arma::ivec& startSamplingPeriod,
                           arma::ivec& endSamplingPeriod,
                           arma::field<arma::imat>& TestField, 
                           arma::field<arma::ivec>& TestTimes,
                           arma::vec& hp_theta,
                           arma::vec& hp_rho,
                           arma::mat Sigma2){
  
  int numTests = thetas.n_elem;
  arma::vec out = arma::zeros<arma::vec>(2*numTests);

  arma::vec cur = arma::zeros<arma::vec>(2*numTests);
  for (int iTest=0; iTest<numTests; iTest++) {
    cur[iTest] = thetas[iTest];
    cur[iTest+numTests] = rhos[iTest];
  }
  
  arma::vec curLogit = logit(cur);
  
  // multivariate normal proposal
  arma::vec canLogit = multrnorm(curLogit, Sigma2);
  
  arma::vec thetas_canLogit = arma::zeros<arma::vec>(numTests);
  arma::vec rhos_canLogit = arma::zeros<arma::vec>(numTests);
  for (int iTest=0; iTest<numTests; iTest++) {
    thetas_canLogit[iTest] = canLogit[iTest];
    rhos_canLogit[iTest] = canLogit[iTest+numTests];
  }

  arma::vec thetas_can = logistic(thetas_canLogit);
  arma::vec rhos_can = logistic(rhos_canLogit);
  
  double logPost_can = logPostThetasRhos(thetas_can, 
                                         rhos_can, 
                                         X, 
                                         startSamplingPeriod,
                                         endSamplingPeriod,
                                         TestField,
                                         TestTimes,
                                         hp_theta, 
                                         hp_rho);
  
  double logPost_cur = logPostThetasRhos(thetas, 
                                         rhos, 
                                         X, 
                                         startSamplingPeriod,
                                         endSamplingPeriod,
                                         TestField,
                                         TestTimes,
                                         hp_theta, 
                                         hp_rho);
  
  double logPostDiff = logPost_can - logPost_cur;

  double prob = exp(logPostDiff);
  
  arma::vec toCompare = {1.0, prob};
  double alpha = min(toCompare);
  double u = runif(1,0,1)[0];

  if (u < alpha){
    out = logistic(canLogit);
  } else {
    out = cur;
  }

  return(out);
}
