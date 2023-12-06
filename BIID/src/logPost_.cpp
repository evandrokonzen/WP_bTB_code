#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
double logPost_(arma::vec& logPars, int G, 
                arma::imat& X, 
                arma::imat& totalNumInfec,
                arma::imat& SocGroup,
                arma::imat& totalmPerGroup,
                arma::ivec& birthTimes,
                arma::ivec& startSamplingPeriod,
                arma::ivec& lastObsAliveTimes, 
                arma::imat& capturesAfterMonit,
                arma::imat& ageMat, 
                arma::vec& hp_lambda,
                arma::vec& hp_beta,
                arma::vec& hp_q,
                arma::vec& hp_tau,
                arma::vec& hp_a2,
                arma::vec& hp_b2,
                arma::vec& hp_c1, 
                int k, 
                double K) {
  
  int m = X.n_rows;
  
  double a;
  
  double lambda = exp(logPars[G]);
  arma::vec alpha_js = arma::zeros<arma::vec>(G);
  for(int g=0; g<G; g++){
    alpha_js[g] = exp(logPars[g])*lambda;
  }
  double b = exp(logPars[G+1]);
  double q = logisticD(logPars[G+2]);
  double ql = logPars[G+2];
  double tau = exp(logPars[G+3]);
  double a2 = exp(logPars[G+4]);
  double b2 = exp(logPars[G+5]);
  double c1 = exp(logPars[G+6]);
  
  double logpost;
  double loglik = 0.0;
  
  for(int i=0; i<m; i++){
    
    
    int mint_i = startSamplingPeriod[i];
    
    // int mint_i;
    // if(birthTimes[i]>1L){
    //   mint_i = birthTimes[i];
    // }else{
    //   mint_i = 1L;
    // }
    
    
    if(birthTimes[i] < startSamplingPeriod[i]){
      loglik += logS(ageMat(i, startSamplingPeriod[i] - birthTimes[i]), a2, b2, c1);
    }
    
    for(int j=mint_i; j<lastObsAliveTimes[i]; j++){
      
      int g = SocGroup(i, j-1);
      // if(g==0L){
      //   Rcout << "g at logpost = " << g << std::endl;
      //   Rcpp::stop("found g=0.");
      // }
      
      double age_ij = (double)(ageMat(i,j));
      double log_pti = TrProbSurvive_(age_ij, a2, b2, c1, true);
      int z_t_1 = X(i,j-1);
      int z_t = X(i,j);
      
      if( ((z_t_1==0L) || (z_t_1==1L) || (z_t_1==3L)) && (z_t==9L) ){
        double log_qti = TrProbDeath_(age_ij, a2, b2, c1, true);
        loglik += log_qti;
      }else if((z_t_1==0L) && (z_t==0L)){
        double inf_mgt = totalNumInfec(g-1, j-1)/(pow((double)(totalmPerGroup(g-1, j-1))/K, q));
        a = alpha_js[g-1];
        loglik += log_pti - a - b*inf_mgt;
      }else if((z_t_1==0L) && (z_t==3L)){
        double inf_mgt = totalNumInfec(g-1, j-1)/(pow((double)(totalmPerGroup(g-1, j-1))/K, q));
        a = alpha_js[g-1];
        // loglik += log_pti + log(1-exp(-a-b*inf_mgt));
        loglik += log_pti + Rf_log1mexp(a+b*inf_mgt);
      }else if((z_t_1==3L) && (z_t==3L)){
        // loglik += log_pti - 1/tau;
        loglik += log_pti + log(1 - ErlangCDF(1, k, tau/((double)k)));
      }else if((z_t_1==3L) && (z_t==1L)){
        // loglik += log_pti + log(1-exp(-1/tau));
        // loglik += log_pti + Rf_log1mexp(1/tau);
        loglik += log_pti + log(ErlangCDF(1, k, tau/((double)k)));
      }else if((z_t_1==1L) && (z_t==1L)){
        loglik += log_pti;
      }
      
    }
  }
  
  
  // add correction term for the captures occurring after the monitoring period
  int i;
  int lastCaptTime;
  int numRows = capturesAfterMonit.n_rows;
  for(int ir=0; ir<numRows; ir++){
    
    i = capturesAfterMonit(ir, 0)-1L;
    lastCaptTime = capturesAfterMonit(ir, 1);
    
    for(int j=lastObsAliveTimes[i]; j<lastCaptTime; j++){
      double age_ij = (double)(ageMat(i,j));
      double log_pti = TrProbSurvive_(age_ij, a2, b2, c1, true); 
      loglik += log_pti;
    }
    
  }
    
  
  
  
  // double a_prior = 0.0;
  // for(int g=0; g<G; g++){
  //   a = exp(logPars[g]);
  //   a_prior += Rf_dgamma(a, 1.0, 1.0, 1L) + log(a);
  // }
  // double lambda_prior = Rf_dgamma(lambda, hp_lambda[0], 1.0/hp_lambda[1], 1L) + log(lambda);
  
  double a_prior = 0.0;
  for(int g=0; g<G; g++){
    a_prior += - exp(logPars[g]) + logPars[g];
  }
  double lambda_prior = - hp_lambda[1]*lambda  + log(lambda);
  
  double b_prior = Rf_dgamma(b, hp_beta[0], 1.0/hp_beta[1], 1L) + log(b);
  double q_prior = hp_q[0]*ql - (hp_q[0] + hp_q[1])*log(1 + exp(ql));
  double tau_prior = Rf_dgamma(tau, hp_tau[0], 1.0/hp_tau[1], 1L) + log(tau);
  double a2_prior = Rf_dgamma(a2, hp_a2[0], 1.0/hp_a2[1], 1L) + log(a2);
  double b2_prior = Rf_dgamma(b2, hp_b2[0], 1.0/hp_b2[1], 1L) + log(b2);
  double c1_prior = Rf_dgamma(c1, hp_c1[0], 1.0/hp_c1[1], 1L) + log(c1);
  
  double logprior = a_prior + lambda_prior + b_prior + q_prior + tau_prior +
    a2_prior + b2_prior + c1_prior;

  
  logpost = loglik + logprior;
  
  return(logpost);
}
