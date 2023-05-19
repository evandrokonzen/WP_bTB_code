#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::rowvec normTransProbRest(arma::rowvec& logProbs) {  
  
  int n = logProbs.n_elem;
  double B = max(logProbs);
  double lse = B + log(sum(exp(logProbs - B)));
  arma::rowvec out(n);
  for (int j=0; j<n; j++) {
    out[j] = exp(logProbs[j] - lse);
  }

  return(out);
}
