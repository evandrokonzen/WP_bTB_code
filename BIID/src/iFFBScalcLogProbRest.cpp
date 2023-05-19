#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
void iFFBScalcLogProbRest(int i,
                          int ttt,
                          arma::cube& logProbRest,
                          arma::imat& X,
                          arma::imat& SocGroup,
                          arma::mat& LogProbDyingMat, 
                          arma::mat& LogProbSurvMat,
                          arma::mat& logProbStoSgivenSorE, 
                          arma::mat& logProbStoEgivenSorE, 
                          arma::mat& logProbStoSgivenI, 
                          arma::mat& logProbStoEgivenI, 
                          arma::mat& logProbStoSgivenD, 
                          arma::mat& logProbStoEgivenD, 
                          double& logProbEtoE, 
                          double& logProbEtoI) {

  // bool eval = true;
  int g = SocGroup(i, ttt);
  // if(g==0L){
  //   eval = false;
  // }
  // if(X(i, ttt)==9L){
  //   eval = false;
  // }
  
  // if(eval){

    int state_t = X(i, ttt);
    int state_t1 = X(i, ttt+1);

    if((state_t==0) && (state_t1==0)){
      logProbRest(ttt,0,i) = LogProbSurvMat(i, ttt+1) + logProbStoSgivenSorE(g-1L, ttt);
      logProbRest(ttt,1,i) = LogProbSurvMat(i, ttt+1) + logProbStoSgivenSorE(g-1L, ttt);
      logProbRest(ttt,2,i) = LogProbSurvMat(i, ttt+1) + logProbStoSgivenI(g-1L, ttt);
      logProbRest(ttt,3,i) = LogProbSurvMat(i, ttt+1) + logProbStoSgivenD(g-1L, ttt);
    }else if((state_t==0) && (state_t1==3)){
      logProbRest(ttt,0,i) = LogProbSurvMat(i, ttt+1) + logProbStoEgivenSorE(g-1L, ttt);
      logProbRest(ttt,1,i) = LogProbSurvMat(i, ttt+1) + logProbStoEgivenSorE(g-1L, ttt);
      logProbRest(ttt,2,i) = LogProbSurvMat(i, ttt+1) + logProbStoEgivenI(g-1L, ttt);
      logProbRest(ttt,3,i) = LogProbSurvMat(i, ttt+1) + logProbStoEgivenD(g-1L, ttt);
    }else if((state_t==3) && (state_t1==3)){
      logProbRest(ttt,0,i) = LogProbSurvMat(i, ttt+1) + logProbEtoE;
      logProbRest(ttt,1,i) = LogProbSurvMat(i, ttt+1) + logProbEtoE;
      logProbRest(ttt,2,i) = LogProbSurvMat(i, ttt+1) + logProbEtoE;
      logProbRest(ttt,3,i) = LogProbSurvMat(i, ttt+1) + logProbEtoE;
    }else if((state_t==3) && (state_t1==1)){
      logProbRest(ttt,0,i) = LogProbSurvMat(i, ttt+1) + logProbEtoI;
      logProbRest(ttt,1,i) = LogProbSurvMat(i, ttt+1) + logProbEtoI;
      logProbRest(ttt,2,i) = LogProbSurvMat(i, ttt+1) + logProbEtoI;
      logProbRest(ttt,3,i) = LogProbSurvMat(i, ttt+1) + logProbEtoI;
    }else if((state_t==1) && (state_t1==1)){
      logProbRest(ttt,0,i) = LogProbSurvMat(i, ttt+1);
      logProbRest(ttt,1,i) = LogProbSurvMat(i, ttt+1);
      logProbRest(ttt,2,i) = LogProbSurvMat(i, ttt+1);
      logProbRest(ttt,3,i) = LogProbSurvMat(i, ttt+1);
    }else if((state_t==0) && (state_t1==9)){
      logProbRest(ttt,0,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,1,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,2,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,3,i) = LogProbDyingMat(i, ttt+1);
    }else if((state_t==1) && (state_t1==9)){
      logProbRest(ttt,0,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,1,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,2,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,3,i) = LogProbDyingMat(i, ttt+1);
    }else if((state_t==3) && (state_t1==9)){
      logProbRest(ttt,0,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,1,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,2,i) = LogProbDyingMat(i, ttt+1);
      logProbRest(ttt,3,i) = LogProbDyingMat(i, ttt+1);
    }else if((state_t==3) && (state_t1==0)){
      Rcpp::Rcout << "Some E->S transition was found. " << 
        "This is not allowed in SEI model." << 
          " Check individual id = " << i+1L << std::endl;
      Rcpp::stop("Algorithm stopped.");
    }else if((state_t==1) && (state_t1==3)){
      Rcpp::Rcout << "Some I->E transition was found. " << 
        "This is not allowed in SEI model." << 
          " Check individual id = " << i+1L << std::endl;
      Rcpp::stop("Algorithm stopped.");
    }else if((state_t==1) && (state_t1==0)){
      Rcpp::Rcout << "Some I->S transition was found. " << 
        "This is not allowed in SEI model." << 
          " Check individual id = " << i+1L << std::endl;
      Rcpp::stop("Algorithm stopped.");
    }

  // }

}
