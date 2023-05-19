#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec ivecMinus1(arma::ivec v){
  
  arma::uvec out = arma::conv_to<arma::uvec>::from(v-1L);
  
  return(out);
}


// [[Rcpp::export]]
arma::uvec vecSeq(int numTests){
  
  arma::uvec out = arma::linspace<arma::uvec>(3L, 2L+numTests, numTests);
  
  return(out);
}

