#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::ivec MakeSeasonVec_(int numSeasons, int seasonStart, unsigned int maxt){

  arma::uvec seasonsVec = arma::linspace<arma::uvec>(0,(numSeasons-1), numSeasons);
  arma::uvec rows = arma::find(seasonsVec == seasonStart-1L);
  if(rows.n_elem==0){
    Rcpp::stop("seasonStart must be an integer from {1, ..., numSeasons}.");
  }
  
  arma::ivec seasonVec = arma::ones<arma::ivec>(maxt);
  seasonVec[0] = seasonStart;
  for (unsigned int tt=1; tt<maxt; tt++) {
    if(seasonVec[tt-1]<numSeasons){
      seasonVec[tt] = seasonVec[tt-1] + 1L;
    }
  }

  return(seasonVec);
}
