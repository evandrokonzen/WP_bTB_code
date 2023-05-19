#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
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
            arma::ivec& endSamplingPeriod) {
  
  int out;
  
  // range where changes in xi modify likelihood
  int xiMin;
  int xiMax;
  if(can<cur){
    xiMin = can;
    xiMax = cur;
  }else{
    xiMin = cur;
    xiMax = can;
  }
  
  double logpostDiff = logPostXi(xiMin, xiMax, can, hp_xi, TestFieldProposal, TestTimes, 
                                 thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod) - 
                        logPostXi(xiMin, xiMax, cur, hp_xi, TestField, TestTimes,  
                                  thetas, rhos, phis, X, startSamplingPeriod, endSamplingPeriod);

  // double prob = exp(logpostDiff + sum(can) - sum(curLogPars));
  double prob = exp(logpostDiff);
  
  double alpha = std::min(1.0, prob);
  double u = runif(1,0,1)[0];

  if (u < alpha){
    out = can;
    TestField = TestFieldProposal;
  } else {
    out = cur;
    TestFieldProposal = TestField;
  }
  

  return(out);
}
