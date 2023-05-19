#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec multrnorm(arma::vec mu, arma::mat Sigma) {

  Environment pkg = Environment::namespace_env("mgcv");
  Function f = pkg["rmvn"];
  
  NumericVector mu_ = wrap(mu.t());
  NumericMatrix Sigma_ = wrap(Sigma);
    
  NumericVector out_ = f(Named("n")=1L, Named("mu")=mu_, Named("V")=Sigma_);

  arma::vec out = as<arma::vec>(out_);
  
  return out;
}


// [[Rcpp::export]]
double logdmultrnorm(arma::vec x, arma::vec mu, arma::mat Sigma) {
  
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["dmvnorm"];

  NumericVector x_ = wrap(x.t());
  NumericVector mu_ = wrap(mu.t());
  NumericMatrix Sigma_ = wrap(Sigma);
  
  NumericVector out_ = f(Named("x")=x_, Named("mean")=mu_, Named("sigma")=Sigma_, 
                         Named("log")=true);
  
  return out_[0];
}

// [[Rcpp::export]]
arma::vec randu(int n) {
  
  arma::vec out = arma::zeros<arma::vec>(n);
  out.randu();
  
  return out;
}
