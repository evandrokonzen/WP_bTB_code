#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::imat LocateIndiv(const arma::imat& TestMat, arma::ivec& birthTimes){
  
  int maxt = max(TestMat.col(0));
  int m = max(TestMat.col(1));
  
  arma::imat SocGroup = arma::zeros<arma::imat>(m, maxt);

  arma::ivec id_i = TestMat.col(1);
  for (int i=0; i<m; i++) {
    arma::uvec which_i = arma::find(id_i == i+1L);
    arma::imat Tests_i = TestMat.rows(which_i);
    arma::ivec times_i = Tests_i.col(0);     //
    arma::ivec groups_i = Tests_i.col(2);    //
    
    int tt0 = std::max(0L, birthTimes[i]-1L);
    
    int firstcapttime = min(times_i);
    arma::uvec which_row = arma::find(times_i == firstcapttime);
    int g = groups_i[which_row[0]];  // first group it belongs to 
    // which_row[0] because we assume the group is the same within same quarter

    for (int tt=tt0; tt<maxt; tt++) {
      // check if moved to another group
      arma::uvec tt_capt = arma::find(times_i == tt+1L);
      if(tt_capt.n_elem>0){
        //   //    //   
        // //   checking ifany individual was found in multiple social groups within
        // //   the same quarter
        // if(tt_capt.n_elem>1){
        //   // Rcpp::Rcout << "groups_i.elem(tt_capt): " << 
        //   //   (groups_i.elem(tt_capt)).t() << std::endl;
        //   
        //   if(max(groups_i.elem(tt_capt))!=min(groups_i.elem(tt_capt))){
        //     Rcpp::Rcout << "id: " << i+1 << std::endl;
        //     Rcpp::Rcout << "t: " << tt+1L << std::endl;
        //     Rcpp::Rcout << "groups_i: " << groups_i << std::endl;
        //   }
        // }
        //   //    //    
        int newGroup = groups_i[tt_capt[0]];
        if(newGroup!=g){
          g = newGroup;
        }
      }
      SocGroup(i, tt) = g;
    }

  }

  return(SocGroup);
}
  