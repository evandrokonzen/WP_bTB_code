#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
double logPostXi(int xiMin,
                 int xiMax,
                 double xi,
                 arma::vec& hp_xi,
                 arma::field<arma::imat>& TestField_, 
                 arma::field<arma::ivec>& TestTimes,
                 arma::vec& thetas,
                 arma::vec& rhos,
                 arma::vec& phis,
                 arma::imat& X,
                 arma::ivec& startSamplingPeriod,
                 arma::ivec& endSamplingPeriod) {
  
  // int maxt = X.n_cols;
  int m = X.n_rows;
  
  int numBrockTests = 2L;
  arma::uvec idxTests = arma::linspace<arma::uvec>(0,numBrockTests-1L,numBrockTests);
  
  double logLik = 0.0;
  for (int jj=0; jj<m; jj++) {
    
    int id = jj+1L;
    // int birthTime = birthTimes[jj];
    arma::imat TestMat_i = TestField_(jj);
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
     
     if((tt+t0+1L >= xiMin)&&(tt+t0+1L < xiMax)){
     
      arma::uvec rows = arma::find(TestTimes_i - t0 == tt+1L);
  
      if(rows.n_elem>0){
        
        arma::imat TestMat_i_tt_allTests = TestMat_i.rows(rows);
        arma::imat TestMat_i_tt = TestMat_i_tt_allTests.cols(idxTests); // Brock test columns
  
        for (unsigned int ir=0; ir<rows.n_elem; ir++) {
          arma::ivec Tests_ir = (TestMat_i_tt.row(ir)).t();
          arma::uvec idx = idxTests.elem(arma::find((Tests_ir==0L) || (Tests_ir==1L)));
          
          for (unsigned int ic=0; ic<idx.n_elem; ic++) {
            
              int i = idx[ic];
              
              if(X(id-1,tt+t0)==0L){
                logLik += log(pow(1-phis[i], TestMat_i_tt(ir,i))*
                  pow(phis[i], (1-TestMat_i_tt(ir,i))));
              }else if(X(id-1,tt+t0)==3L){
                logLik += log(pow(thetas[i]*rhos[i], TestMat_i_tt(ir,i))*
                  pow(1-thetas[i]*rhos[i], (1-TestMat_i_tt(ir,i))));
              }else if(X(id-1,tt+t0)==1L){
                logLik += log(pow(thetas[i], TestMat_i_tt(ir,i))*
                  pow(1-thetas[i], (1-TestMat_i_tt(ir,i))));
              }
          
          }
          
        }
  
      }
    } // end if
    }
  
  }
  
  double xiLogPrior = Rf_dnorm4(xi, hp_xi[0], hp_xi[1], 1L);

  double logPost = logLik + xiLogPrior;
  
  return(logPost);
}
