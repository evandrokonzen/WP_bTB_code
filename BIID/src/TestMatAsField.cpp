#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::field<arma::imat> TestMatAsField(const arma::imat& TestMat, int m){
  
  int numTests = TestMat.n_cols - 3L;
  arma::uvec cols = arma::linspace<arma::uvec>(3L, 2L+numTests, numTests);
  
  arma::field<arma::imat> F(m);
  
  arma::ivec id_i = TestMat.col(1);
  for (int i=0; i<m; i++) {
    arma::uvec which_i = arma::find(id_i == i+1L);
    arma::imat Tests_i = TestMat.rows(which_i);
    F(i) = Tests_i.cols(cols);
  }
  
  return(F);
}

// [[Rcpp::export]]
void TestMatAsFieldProposal(arma::field<arma::imat>& TestFieldProposal, 
                            const arma::field<arma::imat>& TestField,
                            const arma::field<arma::ivec>& TestTimes,
                            int xi, int xiCan, int m){

  int numCapt;
  int t;    // t in natural values t=1,2,...,maxt as in TestMat
  int brock1;
  for (int i=0; i<m; i++) {
    
    arma::ivec TestTimes_i = TestTimes(i);
    arma::imat Tests_i = TestFieldProposal(i);
    numCapt = Tests_i.n_rows;

    // Rcpp::Rcout << "i = " << i << std::endl;
    
    // if(i==1000-1){
    //   Rcpp::Rcout << "i = " << i << std::endl;
    //   Rcpp::Rcout << "TestField.row( i ) = " << TestField.row( i ) << std::endl;
    //   Rcpp::Rcout << "TestFieldProposal.row( i ) = " << TestFieldProposal.row( i ) << std::endl;
    //   Rcpp::Rcout << "-------   \n " << std::endl;
    // }
    
    
    
    
    if(xiCan<xi){  // proposing an earlier changepoint
      for (int irow=0; irow<numCapt; irow++) {
        t = TestTimes_i[irow]; 
        if( (t>=xiCan) && (t<xi) ){
          brock1 = Tests_i(irow, 0);
          Tests_i(irow, 0) = Tests_i(irow, 1);
          Tests_i(irow, 1) = brock1;
        }
      }
    }else{  // proposing a later changepoint
      for (int irow=0; irow<numCapt; irow++) {
        t = TestTimes_i[irow]; 
        if( (t>=xi) && (t<xiCan) ){
          brock1 = Tests_i(irow, 0);
          Tests_i(irow, 0) = Tests_i(irow, 1);
          Tests_i(irow, 1) = brock1;
        }
      }
    }
    
    TestFieldProposal(i) = Tests_i;
      
      // if(i==1000-1){
      //   Rcpp::Rcout << "After correction: " << std::endl;
      //   Rcpp::Rcout << "TestField.row( i ) = " << TestField.row( i ) << std::endl;
      //   Rcpp::Rcout << "TestFieldProposal.row( i ) = " << TestFieldProposal.row( i ) << std::endl;
      //   Rcpp::Rcout << "-------   \n " << std::endl;
      // }
      
  }
  
}
