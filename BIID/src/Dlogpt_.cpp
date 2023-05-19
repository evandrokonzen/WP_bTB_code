#include <Rcpp.h>

// [[Rcpp::export]]
double Dlogpt_a2(double age, double a2, double b2){
  // diffExpsLateLife = exp(b2*(age-1)) - exp(b2*age)
  double y1 = b2*(age-1);
  double y2 = b2*age;
  double diffExpsLateLife = -exp(y1 + log(exp(y2-y1)-1));
  return((a2/b2)*diffExpsLateLife);
}

// [[Rcpp::export]]
double Dlogpt_b2(double age, double a2, double b2){
  
  // diffExpsLateLife = exp(b2*(age-1)) - exp(b2*age)
  double y1 = b2*(age-1);
  double y2 = b2*age;
  double diffExpsLateLife = -exp(y1 + log(exp(y2-y1)-1));
  
  // diffExpsLateLife2 = ((age-1)*exp(b2*(age-1)) - age*exp(b2*age))
  double x1 = b2*(age-1) + log(age-1);
  double x2 = b2*age + log(age);
  double diffExpsLateLife2 = -exp(x1 + log(exp(x2-x1)-1));
  if(age==1.0){
    diffExpsLateLife2 = -age*exp(b2*age);
    // diffExpsLateLife2 = ((age-1)*exp(b2*(age-1)) - age*exp(b2*age));
  }
  return((-a2/b2) * diffExpsLateLife + a2*diffExpsLateLife2);
}

// [[Rcpp::export]]
double Dlogpt_c1(double c1){
  return(-c1);
}
