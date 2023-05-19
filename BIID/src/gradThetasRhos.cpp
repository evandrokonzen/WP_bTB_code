#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::vec gradThetasRhos(arma::vec& thetas,
                         arma::vec& rhos,
                         arma::imat& X,
                         arma::ivec& startSamplingPeriod,
                         arma::ivec& endSamplingPeriod,
                         arma::field<arma::imat>& TestField, 
                         arma::field<arma::ivec>& TestTimes, 
                         arma::vec& hp_theta,
                         arma::vec& hp_rho){
  
  // int maxt = X.n_cols;
  int m = X.n_rows;

  int numTests = TestField(0).n_cols;
  arma::uvec idxTests = arma::linspace<arma::uvec>(0,numTests-1L,numTests);

  arma::vec derivloglik = arma::zeros<arma::vec>(2*numTests);
  
  for (int jj=0; jj<m; jj++) {
  
    int id = jj+1L;
    // int birthTime = birthTimes[jj];
    arma::imat TestMat_i = TestField(jj);
    arma::ivec TestTimes_i = TestTimes(jj);
  
    int t0 = startSamplingPeriod[jj] - 1L;
    int maxt_i = endSamplingPeriod[jj] - t0;
  
    // int t0;
    // int maxt_i;
    // if(birthTime<1L){
    //   t0 = 0L;
    //   // maxt_i = maxt;
    //   maxt_i = endSamplingPeriod[jj];
    // }else{
    //   t0 = birthTime-1L;
    //   // maxt_i = maxt - t0;
    //   maxt_i = endSamplingPeriod[jj] - t0;
    // }

    for (int tt=0; tt<maxt_i; tt++) {
     
      arma::uvec rows = arma::find(TestTimes_i - t0 == tt+1L);
  
      if(rows.n_elem>0){
        
        arma::imat TestMat_i_tt = TestMat_i.rows(rows);
  
        for (unsigned int ir=0; ir<rows.n_elem; ir++) {
          arma::ivec Tests_ir = (TestMat_i_tt.row(ir)).t();
          arma::uvec idx = idxTests.elem(arma::find((Tests_ir==0L) || (Tests_ir==1L)));
          
          for (unsigned int ic=0; ic<idx.n_elem; ic++) {
            
              int i = idx[ic];
              if(X(id-1,tt+t0)==3L){
                double expThetaTilde = exp(logitD(thetas[i]));
                double expRhoTilde = exp(logitD(rhos[i]));
                
                // derivatives wrt theta
                derivloglik[i] += TestMat_i_tt(ir,i)*(1 - thetas[i]) + 
                             (1 - TestMat_i_tt(ir,i))*(
                    expThetaTilde/(1+expThetaTilde+expRhoTilde) - thetas[i]);

                // derivatives wrt rho
                derivloglik[i+numTests] += TestMat_i_tt(ir,i)*(1 - rhos[i]) + 
                                      (1 - TestMat_i_tt(ir,i))*(
                    expRhoTilde/(1+expThetaTilde+expRhoTilde) - rhos[i]);
              }else if(X(id-1,tt+t0)==1L){
                // derivatives wrt theta
                derivloglik[i] += TestMat_i_tt(ir,i)*(1 - thetas[i]) -  
                  (1 - TestMat_i_tt(ir,i))*thetas[i];
              }
          
          }
          
        }
  
      }
      
    }

  }


  arma::vec derivLogPriorWithJac = arma::zeros<arma::vec>(2*numTests);
  for (int iTest=0; iTest<numTests; iTest++) {
    derivLogPriorWithJac[iTest] = hp_theta[0] - (hp_theta[0]+hp_theta[1])*thetas[iTest];
    derivLogPriorWithJac[iTest+numTests] = hp_rho[0] - (hp_rho[0]+hp_rho[1])*rhos[iTest];
  }
  
  arma::vec grad = derivloglik + derivLogPriorWithJac;

  return(grad);
}
