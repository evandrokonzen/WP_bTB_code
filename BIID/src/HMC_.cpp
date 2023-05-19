#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::vec HMC_(arma::vec curLogPars, int G, 
               arma::imat& X, 
               arma::imat& totalNumInfec,
               arma::imat& SocGroup,
               arma::imat& totalmPerGroup,
               arma::ivec& birthTimes, 
               arma::ivec& startSamplingPeriod,
               arma::ivec& lastObsAliveTimes,
               arma::imat& capturesAfterMonit,
               arma::imat& ageMat,
               double epsilon, double L, 
               arma::vec& hp_lambda,
               arma::vec& hp_beta,
               arma::vec& hp_q,
               arma::vec& hp_tau,
               arma::vec& hp_a2,
               arma::vec& hp_b2,
               arma::vec& hp_c1,
               int k, double K) {
  
  arma::vec out = arma::zeros<arma::vec>(curLogPars.n_elem);
  
  arma::vec q = curLogPars;

  arma::vec p = as<arma::vec>(rnorm(q.n_elem, 0, 1));

  arma::vec curp = p;

  p = p + epsilon * grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                          birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                          hp_lambda, hp_beta, hp_q, hp_tau,
                          hp_a2, hp_b2, hp_c1, k, K)/2;
  
  // Rcout << "p first adj = " << p.t() << std::endl;

  int intL = ceil(runif(1,0,1)[0]*L);
  // int intL = L;
  
  for(int i=0; i<intL-1L; i++){
    
    // Rcout << "i = " << i << " out of L=" << intL-1L-1L << std::endl;

    q += epsilon*p;
    
    p += epsilon * grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                         birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                         hp_lambda, hp_beta, hp_q, hp_tau,
                         hp_a2, hp_b2, hp_c1, k, K);

  }
  q = q + epsilon*p;
  p += epsilon * grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                       birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                       hp_lambda, hp_beta, hp_q, hp_tau,
                       hp_a2, hp_b2, hp_c1, k, K)/2;
  p = -p;

  double ProposedH = logPost_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                              birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                              hp_lambda, hp_beta, hp_q, hp_tau,
                              hp_a2, hp_b2, hp_c1, k, K) - 0.5*arma::dot(p,p);
  double CurrentH = logPost_(curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                             birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                             hp_lambda, hp_beta, hp_q, hp_tau,
                             hp_a2, hp_b2, hp_c1, k, K) - 0.5*arma::dot(curp,curp);
  
  // double prob = exp(ProposedH + sum(q) - CurrentH - sum(curLogPars));
  double prob = exp(ProposedH - CurrentH);
  
  double alpha = std::min(1.0, prob);
  double u = runif(1,0,1)[0];
  
  if (u < alpha){
    out = q;
  } else {
    out = curLogPars;
  }
  return(out);
}

// double alpha = std::min(1.0, exp(ap));


// [[Rcpp::export]]
arma::vec HMC_2(arma::vec curLogPars, int G, 
               arma::imat& X, 
               arma::imat& totalNumInfec,
               arma::imat& SocGroup,
               arma::imat& totalmPerGroup,
               arma::ivec& birthTimes, 
               arma::ivec& startSamplingPeriod,
               arma::ivec& lastObsAliveTimes,
               arma::imat& capturesAfterMonit,
               arma::imat& ageMat,
               double epsilon, 
               double epsilonalphas, 
               double epsilonbq, 
               double epsilontau,
               double epsilonc1, 
               int nParsNotGibbs, double L, 
               arma::vec& hp_lambda,
               arma::vec& hp_beta,
               arma::vec& hp_q,
               arma::vec& hp_tau,
               arma::vec& hp_a2,
               arma::vec& hp_b2,
               arma::vec& hp_c1,
               int k, double K) {
  
  
  arma::vec epsilon_ = arma::zeros<arma::vec>(nParsNotGibbs);
  
  for(int j=0; j<nParsNotGibbs; j++){
    if(j<G){
      epsilon_[j] = epsilonalphas;
    }else if((j==G+1) || (j==G+2)){
       epsilon_[j] = epsilonbq;
    }else if(j==G+3){
      epsilon_[j] = epsilontau;
    }else if(j==G+6){
      epsilon_[j] = epsilonc1;
    }else{
      epsilon_[j] = epsilon;
    }
  }
    

  arma::vec out = arma::zeros<arma::vec>(curLogPars.n_elem);
  
  arma::vec q = curLogPars;
  
  arma::vec p = as<arma::vec>(rnorm(q.n_elem, 0, 1));
  
  arma::vec curp = p;
  
  p = p + epsilon_ % grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                          birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                          hp_lambda, hp_beta, hp_q, hp_tau,
                          hp_a2, hp_b2, hp_c1, k, K)/2;
  
  // Rcout << "p first adj = " << p.t() << std::endl;
  
  int intL = ceil(runif(1,0,1)[0]*L);
  // int intL = L;
  
  for(int i=0; i<intL-1L; i++){
    
    // Rcout << "i = " << i << " out of L=" << intL-1L-1L << std::endl;
    
    q += epsilon_ % p;
    
    p += epsilon_ % grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                         birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                         hp_lambda, hp_beta, hp_q, hp_tau,
                         hp_a2, hp_b2, hp_c1, k, K);
    
  }
  q = q + epsilon_ % p;
  p += epsilon_ % grad_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                       birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                       hp_lambda, hp_beta, hp_q, hp_tau,
                       hp_a2, hp_b2, hp_c1, k, K)/2;
  p = -p;
  
  double ProposedH = logPost_(q, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                              birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                              hp_lambda, hp_beta, hp_q, hp_tau,
                              hp_a2, hp_b2, hp_c1, k, K) - 0.5*arma::dot(p,p);
  double CurrentH = logPost_(curLogPars, G, X, totalNumInfec, SocGroup, totalmPerGroup,
                             birthTimes, startSamplingPeriod, lastObsAliveTimes, capturesAfterMonit, ageMat, 
                             hp_lambda, hp_beta, hp_q, hp_tau,
                             hp_a2, hp_b2, hp_c1, k, K) - 0.5*arma::dot(curp,curp);
  
  // double prob = exp(ProposedH + sum(q) - CurrentH - sum(curLogPars));
  double prob = exp(ProposedH - CurrentH);
  
  double alpha = std::min(1.0, prob);
  double u = runif(1,0,1)[0];
  
  if (u < alpha){
    out = q;
  } else {
    out = curLogPars;
  }
  return(out);
}

// double alpha = std::min(1.0, exp(ap));