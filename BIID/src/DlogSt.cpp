#include <Rcpp.h>

// [[Rcpp::export]]
double logS(double age, double a2, double b2, double c1){
  return( - c1*age + (a2/b2)*(1-exp(b2*age)));
}

// [[Rcpp::export]]
double DlogS_a2(double age, double a2, double b2){
  return((a2/b2)*(1-exp(b2*age)));
}

// [[Rcpp::export]]
double DlogS_b2(double age, double a2, double b2){
  return((-a2/b2)*(1 - exp(b2*age)  + exp(b2*age)*b2*age));
}

// [[Rcpp::export]]
double DlogS_c1(double age, double c1){
  return(-c1*age);
}
