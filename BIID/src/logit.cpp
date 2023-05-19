#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
double logitD(double x){
  return(log( x / (1 - x) ));
}

// [[Rcpp::export]]
double logisticD(double x){
  return(exp(x) / (1 + exp(x)));
}

// [[Rcpp::export]]
arma::vec logit(arma::vec& x){
  int n = x.n_elem;
  arma::vec logit = arma::zeros<arma::vec>(n);
  for (int i=0; i<n; i++) {
    logit[i] = log( x[i] / (1 - x[i]) );
  }
  return(logit);
}

// [[Rcpp::export]]
arma::vec logistic(arma::vec& x){
  int n = x.n_elem;
  arma::vec logistic = arma::zeros<arma::vec>(n);
  for (int i=0; i<n; i++) {
    logistic[i] = exp(x[i]) / (1 + exp(x[i]));
  }
  return(logistic);
}

// [[Rcpp::export]]
double sumLogJacobian(arma::vec& xtilde){
  int n = xtilde.n_elem;
  double sumLogsJac = 0.0;
  for (int i=0; i<n; i++) {
    sumLogsJac += xtilde[i] - 2*log(1+exp(xtilde[i]));
  }
  return(sumLogsJac);
}

