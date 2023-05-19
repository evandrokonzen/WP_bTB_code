#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
void ObsProcess_(arma::mat& corrector,
                 int t0,
                 int endTime,
                 int id,
                 const arma::imat& CaptHist, 
                 arma::imat& TestMat_i,
                 arma::ivec& TestTimes_i,
                 arma::vec& etas, 
                 arma::vec& thetas,
                 arma::vec& rhos,
                 arma::vec& phis, 
                 arma::ivec& seasonVec) {
  
  int numTests = TestMat_i.n_cols;
  
  arma::uvec idxTests = arma::linspace<arma::uvec>(0,numTests-1L,numTests);

  // int maxt_i = CaptHist.n_cols - t0;
  int maxt_i = endTime - t0;

  for (int tt=0; tt<maxt_i; tt++) {
   
   double eta = etas[seasonVec[tt+t0] - 1L];
    
    if(CaptHist(id-1L, tt+t0)==0L){
      corrector.row(tt+t0) = {1-eta, 1-eta, 1-eta, 1};
    }else{
      corrector.row(tt+t0) = {eta, eta, eta, 0};
      
      arma::uvec rows = arma::find(TestTimes_i - t0 == tt+1L);

      if(rows.n_elem>0){
        
        arma::imat TestMat_i_tt = TestMat_i.rows(rows);

        double productIfSuscep = 1.0;
        double productIfExposed = 1.0;
        double productIfInfectious = 1.0;
        
        for (unsigned int ir=0; ir<rows.n_elem; ir++) {
          arma::ivec Tests_ir = (TestMat_i_tt.row(ir)).t();
          arma::uvec idx = idxTests.elem(arma::find((Tests_ir==0L) || (Tests_ir==1L)));
          for (unsigned int ic=0; ic<idx.n_elem; ic++) {
            int i = idx[ic];
            productIfSuscep = productIfSuscep*(pow(1-phis[i], TestMat_i_tt(ir,i))*
              pow(phis[i], (1-TestMat_i_tt(ir,i))));
            productIfExposed = productIfExposed*(pow(thetas[i]*rhos[i], TestMat_i_tt(ir,i))*
              pow(1-thetas[i]*rhos[i], (1-TestMat_i_tt(ir,i))));
            productIfInfectious = productIfInfectious*(pow(thetas[i], TestMat_i_tt(ir,i))*
              pow(1-thetas[i], (1-TestMat_i_tt(ir,i))));
          }
        }
        corrector(tt+t0,0) *= productIfSuscep;
        corrector(tt+t0,1) *= productIfExposed;
        corrector(tt+t0,2) *= productIfInfectious;

      }
    }
  }
}
