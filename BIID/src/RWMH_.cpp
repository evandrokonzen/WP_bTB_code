#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
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
                int k, double K) {
  
  arma::vec out(curLogPars.n_elem);

  double logpostDiff = logPost_(can, G, X, totalNumInfec, 
                                SocGroup, totalmPerGroup,
                                birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, 
                                ageMat, 
                                hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K) - 
                        logPost_(curLogPars, G, X, totalNumInfec, 
                                 SocGroup, totalmPerGroup,
                                 birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, 
                                 ageMat, 
                                 hp_lambda, hp_beta, hp_q, hp_tau, hp_a2, hp_b2, hp_c1, k, K);

  // double prob = exp(logpostDiff + sum(can) - sum(curLogPars));
  double prob = exp(logpostDiff);
  
  double alpha = std::min(1.0, prob);
  double u = runif(1,0,1)[0];

  if (u < alpha){
    out = can;
  } else {
    out = curLogPars;
  }
  
  return(out);
}

// double alpha = std::min(1.0, prob);
// 

// // Rcout << "runif: " << runif(1,0,1)[0] << std::endl;
// // RNGScope scope;
// 
// // // using using R function mgcv::rmvn():
// // arma::vec can = multrnorm(wrap(curLogPars.t()), wrap(Sigma));
// // arma::vec can = multrnorm(curLogPars, Sigma);
// // arma::vec mu_ = as<arma::vec>(curLogPars);
// // arma::mat Sigma_ = as<arma::mat>(Sigma);
// // arma::vec yMultRnorm = rnorm(Sigma_.n_cols, 0, 1);
// // arma::vec outMultRnorm = mu_ + arma::chol(Sigma_).t()*yMultRnorm;
// // NumericVector can = wrap(outMultRnorm);
// 
// // Rcout << "can: " << can << std::endl;
// 
// //// alternatively, we can use:
// // NumericVector can = wrap(arma::mvnrnd(curLogPars, Sigma, 1));