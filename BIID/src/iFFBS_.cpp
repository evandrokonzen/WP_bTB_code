#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h> // for RcppArmadillo::sample()
#include "functions.h"
#include <Rcpp/Benchmark/Timer.h>

// [[Rcpp::export]]
void iFFBS_(arma::vec& alpha_js, 
            double b, double q, double tau, int k, double K,
            arma::mat& probDyingMat,
            arma::mat& LogProbDyingMat, 
            arma::mat& LogProbSurvMat,
            arma::cube& logProbRest,
            arma::ivec& nuTimes,
            arma::vec& nuEs,
            arma::vec& nuIs,
            arma::vec& thetas, 
            arma::vec& rhos,
            arma::vec& phis,
            arma::vec& etas, 
            int id,
            int birthTime,
            int startTime,
            int endTime,
            arma::imat& X,
            arma::ivec& seasonVec,
            arma::imat& TestMat_i, 
            arma::ivec& TestTimes_i,
            const arma::imat& CaptHist,
            arma::mat& corrector,
            arma::mat& predProb,
            arma::mat& filtProb,
            arma::mat& logTransProbRest,
            arma::imat& numInfecMat, 
            arma::imat& SocGroup,
            arma::imat& mPerGroup,
            arma::ivec& idVecAll, 
            arma::mat& logProbStoSgivenSorE, 
            arma::mat& logProbStoEgivenSorE, 
            arma::mat& logProbStoSgivenI, 
            arma::mat& logProbStoEgivenI, 
            arma::mat& logProbStoSgivenD, 
            arma::mat& logProbStoEgivenD, 
            double& logProbEtoE, 
            double& logProbEtoI,
            arma::field<arma::ivec>& whichRequireUpdate,
            double& sumLogCorrector){
  
  int m = X.n_rows;
  int maxt = X.n_cols;
  
  int numStates = filtProb.n_cols;
  

  
  int t0 = startTime - 1L;
  int maxt_i = endTime - t0;

  // int t0;
  // if(birthTime<1L){
  //   t0 = 0L;
  // }else{
  //   t0 = birthTime-1L;
  // }
  // int maxt_i = endTime - t0;
  
  
  
  // Rcout << "id: " << id << std::endl;
  // Rcout << "birthTime: " << birthTime << std::endl;
  // Rcout << "endTime: " << endTime << std::endl;
  
  // update corrector
  ObsProcess_(corrector, t0, endTime, id, CaptHist, TestMat_i, TestTimes_i, 
              etas, thetas, rhos, phis, seasonVec);
  
  
  // Rcout << "corrector: " << std::endl;

  
  int idNext;
  if(id<m){
    idNext = id; // for the next individual (not id-1)
  }else{
    idNext = 0L; // go back to 1st individual when id==m
  }

  // Forward Filtering --------------------------------------
  
  double prDeath;
  
  arma::rowvec unnormFiltProb = arma::zeros<arma::rowvec>(numStates);
  arma::rowvec transProbRest = arma::zeros<arma::rowvec>(numStates);
  
  // t=1
  // if(birthTime>0L){
    // nuE = 0.0;
    // nuI = 0.0;
  // }
  
  double nuE_i = 0.0;
  double nuI_i = 0.0;
  
  if(birthTime < startTime){  // born before monitoring started
    int nuIdx = min(arma::find(nuTimes == startTime)); // min is used to convert to int type
    nuE_i = nuEs[nuIdx];
    nuI_i = nuIs[nuIdx];
  }// otherwise, nuE_i = nuI_i = 0.0

  // The grid of forward sweep starts at t0.
  // t0: either the beginning of the study or the date of birth
  // The individual must be alive at t0.
  // If it was born at or after the beginning of the study, it's assumed to be 
  // susceptible at t0. Otherwise, nuE and nuI are used.
  predProb(t0,0) = 1.0-nuE_i-nuI_i;
  predProb(t0,1) = nuE_i;
  predProb(t0,2) = nuI_i;
  predProb(t0,3) = 0.0;
  
  // prDeath = probDyingMat(id-1, t0);
  // predProb(t0,0) = (1-prDeath)*(1-nu);
  // predProb(t0,1) = (1-prDeath)*nu;
  // predProb(t0,2) = prDeath;
  
  if(t0 < maxt-1){
    arma::rowvec logTransProbRest_row = logTransProbRest.row(t0);
    transProbRest = normTransProbRest(logTransProbRest_row);
    for(int s=0; s<numStates; s++){
      unnormFiltProb[s] = corrector(t0,s) * predProb(t0,s) * transProbRest[s];
    }
  }else{
    for(int s=0; s<numStates; s++){
      unnormFiltProb[s] = corrector(t0,s) * predProb(t0,s);
    }
  }
  filtProb.row(t0) = unnormFiltProb / sum(unnormFiltProb);
  
  
  if(maxt_i>2){
    // t=2,...,T-1
    for(int tt=1; tt<(maxt_i-1); tt++){
      
      int g = SocGroup(id-1, tt-1+t0);

      double prDeath = probDyingMat(id-1, tt+t0);
      
      double p00 = (1-prDeath)*exp(logProbStoSgivenSorE(g-1L, tt-1+t0)); 
      double p01 = (1-prDeath)*exp(logProbStoEgivenSorE(g-1L, tt-1+t0));
      double p11 = (1-prDeath)*exp(logProbEtoE);
      double p12 = (1-prDeath)*exp(logProbEtoI);
      double p22 = (1-prDeath);
      
      predProb(tt+t0,0) = p00*(filtProb(tt-1+t0,0));
      predProb(tt+t0,1) = p01*(filtProb(tt-1+t0,0)) + p11*filtProb(tt-1+t0,1);
      predProb(tt+t0,2) = p12*(filtProb(tt-1+t0,1)) + p22*filtProb(tt-1+t0,2);
      predProb(tt+t0,3) = prDeath*filtProb(tt-1+t0,0) + 
                          prDeath*filtProb(tt-1+t0,1) + 
                          prDeath*filtProb(tt-1+t0,2) +
                          1*filtProb(tt-1+t0,3);
      
      arma::rowvec logTransProbRest_row = logTransProbRest.row(tt+t0);
      transProbRest = normTransProbRest(logTransProbRest_row);
      for(int s=0; s<numStates; s++){
        unnormFiltProb[s] = corrector(tt+t0,s) * predProb(tt+t0,s) * transProbRest[s];
      }
      
      filtProb.row(tt+t0) = unnormFiltProb / sum(unnormFiltProb);

    }
  }
  
  // t=T
  if(maxt_i>=1){
    
    int tt = maxt_i-1;
    
    int g = SocGroup(id-1, tt-1+t0);
    
    prDeath = probDyingMat(id-1, tt+t0);
    
    double p00 = (1-prDeath)*exp(logProbStoSgivenSorE(g-1L, tt-1+t0)); 
    double p01 = (1-prDeath)*exp(logProbStoEgivenSorE(g-1L, tt-1+t0));
    double p11 = (1-prDeath)*exp(logProbEtoE);
    double p12 = (1-prDeath)*exp(logProbEtoI);
    double p22 = (1-prDeath);
    
    predProb(tt+t0,0) = p00*(filtProb(tt-1+t0,0));
    predProb(tt+t0,1) = p01*(filtProb(tt-1+t0,0)) + p11*filtProb(tt-1+t0,1);
    predProb(tt+t0,2) = p12*(filtProb(tt-1+t0,1)) + p22*filtProb(tt-1+t0,2);
    predProb(tt+t0,3) = prDeath*filtProb(tt-1+t0,0) + 
      prDeath*filtProb(tt-1+t0,1) + 
      prDeath*filtProb(tt-1+t0,2) +
      1*filtProb(tt-1+t0,3);

    if(tt+t0 < maxt-1){
      arma::rowvec logTransProbRest_row = logTransProbRest.row(tt+t0);
      transProbRest = normTransProbRest(logTransProbRest_row);
      for(int s=0; s<numStates; s++){
        unnormFiltProb[s] = corrector(tt+t0,s) * predProb(tt+t0,s) * transProbRest[s];
      }
    }else{  // if(tt+t0 == maxt-1){
      for(int s=0; s<numStates; s++){
        unnormFiltProb[s] = corrector(tt+t0,s) * predProb(tt+t0,s);
      }
    }

    filtProb.row(tt+t0) = unnormFiltProb / sum(unnormFiltProb);

  }
  
  
  
  // Rcout << "forward sweep: " << std::endl;
  
  // Backward Sampling --------------------------------------
  
  // arma::uvec tFlagged = arma::zeros<arma::uvec>(maxt-1);
  
  arma::ivec states = {0L, 3L, 1L, 9L};
  // arma::vec probs = filtProb.row(maxt-1).t();
  arma::vec probs = filtProb.row(endTime-1).t();
  
  
  // int g = SocGroup(id-1, maxt-1);
  // int oldStatus = X(id-1, maxt-1);
  int newStatus = RcppArmadillo::sample(states, 1, true, probs)[0];
  

  // // -------------------------------------------------------------------
  // THIS IS NO LONGER NEEDED; mPerGroup now no longer includes the individual 
  // that is being updated --------------------
  // // updating mPerGroup  
  // if( ((oldStatus==0L) || (oldStatus==3L) || (oldStatus==1L)) && (newStatus==9L) ){
  //   // tFlagged[maxt-1] = 1L;
  //   mPerGroup(g-1, maxt-1) -= 1L;
  // }
  // if( (oldStatus==9L) && ((newStatus==0L) || (newStatus==3L) || (newStatus==1L)) ){
  //   // tFlagged[maxt-1] = 1L;
  //   mPerGroup(g-1, maxt-1) += 1L;
  // }
  // // --------------------------------------------------------------------

  // X(id-1, maxt-1) = newStatus;
  X(id-1, endTime-1) = newStatus;
  
  // double a;
  
  // tt will start from maxt_i-2 and finish at 0
  if(maxt_i>1){
    for(int tt = (maxt_i-1); tt --> 0;){
      
      int g = SocGroup(id-1, tt+t0);
      int mgt = mPerGroup(g-1, tt+t0);
      
      double a = alpha_js[g-1L];
        
      double inf_mgt = numInfecMat(g-1, tt+t0)/(pow((double)(mgt + 1.0)/K, q));  

      double prDeath = probDyingMat(id-1, tt+1+t0);
      
      double p00 = (1-prDeath)*exp(-a - b*inf_mgt);
      double p01 = (1-prDeath)*(1 - exp(-a - b*inf_mgt));
      // double p11 = (1-prDeath)*exp(-1/tau);
      // double p12 = (1-prDeath)*(1 - exp(-1/tau));
      double p11 = (1-prDeath)*(1.0 - ErlangCDF(1, k, tau/((double)k) ));
      double p12 = (1-prDeath)*ErlangCDF(1, k, tau/((double)k) );
      double p22 = (1-prDeath);
      
      double p09 = prDeath;
      double p19 = prDeath;
      double p29 = prDeath;
      
      double probSuscep_t=0.0;
      double probE_t=0.0;
      double probI_t=0.0;
      double probDead_t=0.0;
      
      if(X(id-1,tt+1+t0) == 0L){
        probSuscep_t = (p00*filtProb(tt+t0, 0))/(predProb(tt+1+t0, 0));
        probE_t = 0.0;
        probI_t = 0.0;
        probDead_t = 0.0;
      }else if(X(id-1,tt+1+t0) == 3L){
        probSuscep_t = (p01*filtProb(tt+t0, 0))/(predProb(tt+1+t0, 1));
        probE_t = (p11*filtProb(tt+t0, 1))/(predProb(tt+1+t0, 1));
        probI_t = 0.0;
        probDead_t = 0.0;
      }else if(X(id-1,tt+1+t0) == 1L){
        probSuscep_t = 0.0;
        probE_t = (p12*filtProb(tt+t0, 1))/(predProb(tt+1+t0, 2));
        probI_t = (p22*filtProb(tt+t0, 2))/(predProb(tt+1+t0, 2));
        probDead_t = 0.0;
      }else if(X(id-1,tt+1+t0) == 9L){
        probSuscep_t = (p09*filtProb(tt+t0, 0))/(predProb(tt+1+t0, 3));
        probE_t = (p19*filtProb(tt+t0, 1))/(predProb(tt+1+t0, 3));
        probI_t = (p29*filtProb(tt+t0, 2))/(predProb(tt+1+t0, 3));
        probDead_t = filtProb(tt+t0, 3)/predProb(tt+1+t0, 3);
      }
      
      probs = {probSuscep_t, probE_t, probI_t, probDead_t};
      
      if((tt==0L)&&(birthTime>=startTime)){
        probs = {1.0, 0.0, 0.0, 0.0};
      }

      // int oldStatus = X(id-1, tt+t0);
      int newStatus = RcppArmadillo::sample(states, 1, true, probs)[0];

      // // -------------------------------------------------------------------
      // THIS IS NO LONGER NEEDED; mPerGroup now no longer includes the individual 
      // that is being updated --------------------
      // // updating mPerGroup  
      // if( ((oldStatus==0L) || (oldStatus==3L) || (oldStatus==1L)) && (newStatus==9L) ){
      //   // tFlagged[tt+t0] = 1L;
      //   mPerGroup(g-1, tt+t0) -= 1L;
      // }
      // if( (oldStatus==9L) && ((newStatus==0L) || (newStatus==3L) || (newStatus==1L)) ){
      //   // tFlagged[tt+t0] = 1L;
      //   mPerGroup(g-1, tt+t0) += 1L;
      // }
      // // --------------------------------------------------------------------
      
      
      X(id-1, tt+t0) = newStatus;
      
      // X(id-1, tt+t0) = as<int>(rbinom(1, 1, probs[1]));

    }
  }
  
  // Rcout << "backward probs calculated: " << std::endl;
  
  //   //    //   //    //   //    //   //  
  // calculating log of density of observation process
  
  for (int tt=0; tt<maxt_i; tt++) {
    
    // Rcout << "tt+t0: " << tt+t0 << std::endl;
    // Rcout << "corrector.row(tt+t0): " << corrector.row(tt+t0) << std::endl;
    
    if(X(id-1, tt+t0) == 0L){
      sumLogCorrector += log(corrector(tt+t0, 0L));
    }else if(X(id-1, tt+t0) == 3L){
      sumLogCorrector += log(corrector(tt+t0, 1L));
    }else if(X(id-1, tt+t0) == 1L){
      sumLogCorrector += log(corrector(tt+t0, 2L));
    }else if(X(id-1, tt+t0) == 9L){
      sumLogCorrector += log(corrector(tt+t0, 3L));
    }
    
  }
  //   //    //   //    //   //    //   //    
  
  
  
  
  
  // Rcout << "now will start updating numInfecMat, etc for next ID " << std::endl;
  
  
  // Updating mPerGroup, numInfecMat, logProbStoSgivenSorE, etc for the next individual (id or 0) -----
  
  for(unsigned int tt=0; tt<(maxt-1L); tt++){  // idNext may be included after endTime
  
    int g = SocGroup(id-1, tt);
    int g_idNext = SocGroup(idNext, tt);

    // Rcout << "id: " << id << std::endl;
    // Rcout << "idNext: " << idNext << std::endl;
    // Rcout << "g: " << g << std::endl;
    // Rcout << "g_idNext: " << g_idNext << std::endl;
    // double a = alpha_js[g_idNext-1L];
    // Rcout << "a: " << a << std::endl;
    
    
    
    if((g==g_idNext)&&(g!=0L)){
      
      // // updating/correcting numInfecMat for the next individual
      int infecToAdd = 0L;
      if(X(id-1,tt)==1L){
        infecToAdd += 1L;
      }
      if(X(idNext,tt)==1L){
        infecToAdd -= 1L;
      }

      if((X(id-1,tt)==1L) || (X(idNext,tt)==1L)){
        // tFlagged[tt] = 1L;
        numInfecMat(g_idNext-1L, tt) += infecToAdd;
      }

      // // updating/correcting mPerGroup for the next individual
      int mToAdd = 0L;
      if((X(id-1,tt)==0L) || (X(id-1,tt)==1L)  || (X(id-1,tt)==3L)){
        mToAdd += 1L;
      }
      if((X(idNext,tt)==0L) || (X(idNext,tt)==1L)  || (X(idNext,tt)==3L)){
        mToAdd -= 1L;
      }
      
      if( ((X(id-1,tt)==0L) || (X(id-1,tt)==1L)  || (X(id-1,tt)==3L)) || 
          ((X(idNext,tt)==0L) || (X(idNext,tt)==1L)  || (X(idNext,tt)==3L)) ){
        // tFlagged[tt] = 1L;
        mPerGroup(g_idNext-1L, tt) += mToAdd;
      }
      
      

      
      if(id<m){
        
        double a = alpha_js[g_idNext-1L];
        
        int mgt;
        double inf_mgt;
        mgt = mPerGroup(g_idNext-1L, tt);
        
        inf_mgt = numInfecMat(g_idNext-1L, tt)/(pow((double)(mgt + 1.0)/K, q));
        logProbStoSgivenSorE(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenSorE(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        inf_mgt = (numInfecMat(g_idNext-1L, tt)+1L)/(pow((double)(mgt + 1.0)/K, q));
        logProbStoSgivenI(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenI(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        inf_mgt = numInfecMat(g_idNext-1L, tt)/(pow((double)mgt/K, q));
        logProbStoSgivenD(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenD(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
      }

      
    }else{
      
      if(g!=0L){
        if(X(id-1,tt)==1L){
          numInfecMat(g-1, tt) += 1L;
        }
        if( (X(id-1,tt)==0L) || (X(id-1,tt)==1L)  || (X(id-1,tt)==3L) ){
          mPerGroup(g-1, tt) += 1L;
        }
      }
      
      
      
      if((id<m)&&(g!=0L)){
        
        double a = alpha_js[g-1L];
        
        int mgt;
        double inf_mgt;
        
        mgt = mPerGroup(g-1L, tt); 
        
        inf_mgt = numInfecMat(g-1L, tt)/(pow((double)mgt/K, q));
        logProbStoSgivenSorE(g-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenSorE(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        // inf_mgt = numInfecMat(g-1L, tt)/(pow((double) mgt, q));
        logProbStoSgivenI(g-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenI(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        // inf_mgt = numInfecMat(g-1L, tt)/(pow((double) mgt, q));
        logProbStoSgivenD(g-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenD(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);    
        
        // inf_mgt = numInfecMat(g-1L, tt)/(pow((double) mgt + 1.0, q));
        // logProbStoSgivenSorE(g-1L, tt) = -a -b*inf_mgt;
        // logProbStoEgivenSorE(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        // 
        // inf_mgt = (numInfecMat(g-1L, tt)+1L)/(pow((double) mgt + 1.0, q));
        // logProbStoSgivenI(g-1L, tt) = -a -b*inf_mgt;
        // logProbStoEgivenI(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        // 
        // inf_mgt = numInfecMat(g-1L, tt)/(pow((double) mgt, q));
        // logProbStoSgivenD(g-1L, tt) = -a -b*inf_mgt;
        // logProbStoEgivenD(g-1L, tt) = Rf_log1mexp(a + b*inf_mgt);        

      }
            
      if(g_idNext!=0L){
        if(X(idNext,tt)==1L){
          numInfecMat(g_idNext-1, tt) -= 1L;
        }
        if( (X(idNext,tt)==0L) || (X(idNext,tt)==1L)  || (X(idNext,tt)==3L) ){
          mPerGroup(g_idNext-1, tt) -= 1L;
        }
      }

      if((id<m)&&(g_idNext!=0L)){
        
        double a = alpha_js[g_idNext-1L];
        
        int mgt;
        double inf_mgt;
        mgt = mPerGroup(g_idNext-1L, tt);
        
        inf_mgt = numInfecMat(g_idNext-1L, tt)/(pow((double)(mgt + 1.0)/K, q));
        logProbStoSgivenSorE(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenSorE(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        inf_mgt = (numInfecMat(g_idNext-1L, tt)+1L)/(pow((double)(mgt + 1.0)/K, q));
        logProbStoSgivenI(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenI(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
        inf_mgt = numInfecMat(g_idNext-1L, tt)/(pow((double)mgt/K, q));
        logProbStoSgivenD(g_idNext-1L, tt) = -a -b*inf_mgt;
        logProbStoEgivenD(g_idNext-1L, tt) = Rf_log1mexp(a + b*inf_mgt);
        
      }

      // if((X(id-1,tt)==1L) || (X(idNext,tt)==1L)){
      //   tFlagged[tt] = 1L;
      // }
      
    }
    
  }
  
  // Rcout << "numInfecMat for next ID calculated " << std::endl;
  
  // at id==m, we don't need to update logTransProbRest, because
  // the next individual (id==1) will be updated separately to 
  // consider new rate values
  if(id<m){
    
    int c = (id-1)*(maxt-1);
    int g_1;

    // current individual (id) will be in logProbRest when updating idNext
    for(unsigned int tt=0; tt<maxt-1L; tt++){
    // for(unsigned int tt=0; tt<endTime-1L; tt++){
      
      // update logProbRest(tt,_,id-1);
      // if(X(id-1L, tt)==0L){
        iFFBScalcLogProbRest(id-1L, tt, logProbRest, X, SocGroup, 
                             LogProbDyingMat, LogProbSurvMat, 
                             logProbStoSgivenSorE, logProbStoEgivenSorE, 
                             logProbStoSgivenI, logProbStoEgivenI, 
                             logProbStoSgivenD, logProbStoEgivenD, 
                             logProbEtoE, logProbEtoI);
      // }
      
      

      // adding logProbRest(tt,_,id-1) to logTransProbRest; and 
      // removing logProbRest(tt,s,idNext)
      // (idNext must not be included in logProbRest when updating idNext)
      for(int s=0; s<numStates; s++){
        logTransProbRest(tt, s) += (logProbRest(tt,s,id-1L) - logProbRest(tt,s,idNext));
        logProbRest(tt,s,idNext) = 0.0;
      }
      
      // update logProbRest(tt,_,idNext) at flagged time points for the 
      // remaining m-2 elements
      // This should be required only for individuals that belong to the 
      // same group as id or idNext

      // int pos = (id-1)*(maxt-1) + tt;
      // arma::ivec whichArma = whichRequireUpdate(pos);
      // for(auto & jj : whichArma){
      
      for(auto & jj : whichRequireUpdate(c + tt)){
      
        if(X(jj, tt)==0L){

          for(int s=0; s<numStates; s++){
            logTransProbRest(tt, s) -= logProbRest(tt,s,jj);
          }
          
          // iFFBScalcLogProbRest(jj, tt, logProbRest, X, SocGroup,
          //                      LogProbDyingMat, LogProbSurvMat,
          //                      logProbStoSgivenSorE, logProbStoEgivenSorE, 
          //                      logProbStoSgivenI, logProbStoEgivenI, 
          //                      logProbStoSgivenD, logProbStoEgivenD, 
          //                      logProbEtoE, logProbEtoI);
              
          g_1 = SocGroup(jj, tt)-1L;
          
          if(X(jj, tt+1)==0L){
            logProbRest(tt,0,jj) = LogProbSurvMat(jj, tt+1) + logProbStoSgivenSorE(g_1, tt);
            logProbRest(tt,1,jj) = LogProbSurvMat(jj, tt+1) + logProbStoSgivenSorE(g_1, tt);
            logProbRest(tt,2,jj) = LogProbSurvMat(jj, tt+1) + logProbStoSgivenI(g_1, tt);
            logProbRest(tt,3,jj) = LogProbSurvMat(jj, tt+1) + logProbStoSgivenD(g_1, tt);
          }else if(X(jj, tt+1)==3L){
            logProbRest(tt,0,jj) = LogProbSurvMat(jj, tt+1) + logProbStoEgivenSorE(g_1, tt);
            logProbRest(tt,1,jj) = LogProbSurvMat(jj, tt+1) + logProbStoEgivenSorE(g_1, tt);
            logProbRest(tt,2,jj) = LogProbSurvMat(jj, tt+1) + logProbStoEgivenI(g_1, tt);
            logProbRest(tt,3,jj) = LogProbSurvMat(jj, tt+1) + logProbStoEgivenD(g_1, tt);
          }
              
          for(int s=0; s<numStates; s++){
            logTransProbRest(tt, s) += logProbRest(tt,s,jj);
          }
          

        }
      }

    }
      
  }
  
  // Rcout << "logTransProbRest for next ID calculated: " << std::endl;

}

// Using this:
// for(auto & jj : whichRequireUpdate(c + tt)){
//   if(X(jj, tt)==0L){
// is equivalent to using the loop and if condition below
// 
// for(int jj=0; jj<m; jj++){
// if(((jj!=id-1L)&&(jj!=idNext))&&
//    ((SocGroup(jj, tt)==SocGroup(id-1, tt))||
//     (SocGroup(jj, tt)==SocGroup(idNext, tt)))&&
// (X(jj, tt)==0L)){


//########### 

// List outList(5);
// outList[0] = Xid;
// outList[1] = predProb;
// outList[2] = filtProb;
// return(outList);