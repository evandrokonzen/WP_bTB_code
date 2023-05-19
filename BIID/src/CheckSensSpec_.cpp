#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::imat CheckSensSpec_(int numTests, 
                          arma::field<arma::imat>& TestField, 
                          arma::field<arma::ivec>& TestTimes,
                          arma::imat& X){
  
  int m = X.n_rows;

  arma::imat out = arma::zeros<arma::imat>(4, numTests);

  for (int iTest=0; iTest<numTests; iTest++) {

    int numInfecTested = 0L;
    int numInfecPositives = 0L;
    int numSuscepTested = 0L;
    int numSuscepNegatives = 0L;
    
    for (int i=0; i<m; i++) {  

      arma::imat Tests_i = TestField(i);
      arma::uvec testTimes_i = arma::conv_to<arma::uvec>::from(TestTimes(i)-1L);
      arma::ivec X_i = (X.row(i)).t();
      arma::ivec status = X_i.elem(testTimes_i);

      arma::ivec tests_i = Tests_i.col(iTest);
      
      // exposed and infectious individuals
      arma::uvec which_ExpInfec = arma::find((status == 3L)||(status == 1L));
      arma::ivec tests_i_inf = tests_i.elem(which_ExpInfec);

      int newInfTes = sum((tests_i_inf==0L) || (tests_i_inf==1L));
      int newInfPos = sum(((status == 3L)||(status == 1L)) && (tests_i==1L));

      numInfecTested += newInfTes;
      numInfecPositives += newInfPos;
      
      // susceptible individuals
      arma::uvec which_suscep = arma::find(status == 0L);
      arma::ivec tests_i_suscep = tests_i.elem(which_suscep);
      
      int newSuscepTes = sum((tests_i_suscep==0L) || (tests_i_suscep==1L));
      int newSuscepPos = sum((status==0L) && (tests_i==0L));
      
      numSuscepTested += newSuscepTes;
      numSuscepNegatives += newSuscepPos;
      
    }
    
    out(0, iTest) = numInfecPositives;
    out(1, iTest) = numInfecTested - numInfecPositives;
    out(2, iTest) = numSuscepNegatives;
    out(3, iTest) = numSuscepTested - numSuscepNegatives;
    
  }
  
  return(out);
}
