#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::vec HMC_thetas_rhos(arma::vec& thetas, 
                          arma::vec& rhos, 
                          arma::imat& X,
                          arma::ivec& startSamplingPeriod,
                          arma::ivec& endSamplingPeriod,
                          arma::field<arma::imat>& TestField, 
                          arma::field<arma::ivec>& TestTimes,
                          arma::vec& hp_theta,
                          arma::vec& hp_rho,
                          double epsilon,
                          double L){
  
  int numTests = thetas.n_elem;
  arma::vec out = arma::zeros<arma::vec>(2*numTests);

  // stacking into a unique vector
  arma::vec cur = arma::zeros<arma::vec>(2*numTests);
  for (int iTest=0; iTest<numTests; iTest++) {
    cur[iTest] = thetas[iTest];
    cur[iTest+numTests] = rhos[iTest];
  }

  // multivariate normal proposal
  arma::vec p = as<arma::vec>(rnorm(2*numTests, 0, 1));
  arma::vec curp = p;
  
  arma::vec q_thetas = arma::zeros<arma::vec>(numTests);
  arma::vec q_rhos = arma::zeros<arma::vec>(numTests);
  for (int iTest=0; iTest<numTests; iTest++) {
    q_thetas[iTest] = cur[iTest];
    q_rhos[iTest] = cur[iTest+numTests];
  }

  arma::vec q = logit(cur);
  
  p = p + epsilon * gradThetasRhos(q_thetas, 
                                   q_rhos, 
                                   X, 
                                   startSamplingPeriod,
                                   endSamplingPeriod,
                                   TestField,
                                   TestTimes,
                                   hp_theta, 
                                   hp_rho)/2;
  
  int intL = ceil(runif(1,0,1)[0]*L);                           
  
  for(int i=0; i<intL-1L; i++){

    q = q + epsilon*p;
    
    for (int iTest=0; iTest<numTests; iTest++) {
      q_thetas[iTest] = logisticD(q[iTest]);
      q_rhos[iTest] = logisticD(q[iTest+numTests]);
    }

    p = p + epsilon * gradThetasRhos(q_thetas, 
                                     q_rhos, 
                                     X, 
                                     startSamplingPeriod,
                                     endSamplingPeriod,
                                     TestField,
                                     TestTimes,
                                     hp_theta, 
                                     hp_rho);
    
  }
  q = q + epsilon*p;
  for (int iTest=0; iTest<numTests; iTest++) {
    q_thetas[iTest] = logisticD(q[iTest]);
    q_rhos[iTest] = logisticD(q[iTest+numTests]);
  }
  p = p + epsilon * gradThetasRhos(q_thetas, 
                                   q_rhos, 
                                   X, 
                                   startSamplingPeriod,
                                   endSamplingPeriod,
                                   TestField,
                                   TestTimes,
                                   hp_theta, 
                                   hp_rho)/2;
  p = -p;

  double ProposedH = logPostThetasRhos(q_thetas, q_rhos, X, startSamplingPeriod, endSamplingPeriod,
                                       TestField, TestTimes, 
                                       hp_theta, hp_rho) - 0.5*arma::dot(p,p);
  double CurrentH = logPostThetasRhos(thetas, rhos, X, startSamplingPeriod, endSamplingPeriod,
                                      TestField, TestTimes, 
                                      hp_theta, hp_rho) - 0.5*arma::dot(curp,curp);

  double prob = exp(ProposedH - CurrentH);  

  arma::vec toCompare = {1.0, prob};
  
  double alpha = min(toCompare);
  double u = runif(1,0,1)[0];
  
  if (u < alpha){
    out = logistic(q);
  } else {
    out = cur;
  }

  return(out);
}
