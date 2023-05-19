#include <Rcpp.h>

// [[Rcpp::export]]
double TrProbSurvive_(double age, double a2, double b2, 
                    double c1, bool logar){
  
  // calculating diffExpsLateLife = exp(b2*(age-1)) - exp(b2*age)
  double y1 = b2*(age-1);
  double y2 = b2*age;
  double diffExpsLateLife = -exp(y1 + log(exp(y2-y1)-1));
  double out = - c1 + (a2/b2)*( diffExpsLateLife );
  // double out = - c1 + (a2/b2)*( exp(b2*(age-1)) - exp(b2*age) );
  if(!logar){
    out = exp(out);
  }
  return(out);
}
