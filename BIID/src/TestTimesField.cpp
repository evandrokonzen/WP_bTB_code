#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::field<arma::ivec> TestTimesField(const arma::imat& TestMat, int m){
  
  arma::field<arma::ivec> F(m);
  
  arma::ivec id_i = TestMat.col(1);
  for (int i=0; i<m; i++) {
    arma::uvec which_i = arma::find(id_i == i+1L);
    arma::imat Tests_i = TestMat.rows(which_i);
    F(i) = Tests_i.col(0);
  }
  
  return(F);
}
  