#include <Rcpp.h>

// [[Rcpp::export]]
double fact(int k){
  return(Rf_gammafn((double)(k)+1.0));
}

// [[Rcpp::export]]
double ErlangCDF(int x_, int k, double tau){
  double x = (double)(x_);
  double Q = 1.0;
  if(k>1L){
    for (int n=1L; n<k; n++) {
      Q += (1/fact(n))*pow(x/tau, n);
    }
  }
  double out = 1.0 - exp(-x/tau)*Q;
  return(out);
}

// // example:
// x = 1
// k = 2
// tau = 0.34
// BIID:::ErlangCDF(x_=x, k=k, tau=tau)
// pgamma(q=x, shape=k, scale=tau)
  
// [[Rcpp::export]]
double DerivErlangCDF(int x_, int k, double tau){
  double x = (double)(x_);
  double out = -(1/fact(k-1L))*pow(x/tau, (double)(k))*exp(-x/tau);
  return(out);
}


// // below there's an alternative, equivalent expression
// // [[Rcpp::export]]
// double DerivErlangCDF(int x_, int k, double tau){
//   double x = (double)(x_);
//   double Q1 = 1.0;
//   double Q2 = 0.0;
//   if(k>1L){
//     for (int n=1L; n<k; n++) {
//       Q1 += (1/fact(n))*pow(x/tau, n);
//       Q2 += (1/fact(n-1))*pow(x/tau, n);
//     }
//   }
//   double out = exp(-x/tau)*((-x/tau)*Q1 + Q2);
//   return(out);
// }
