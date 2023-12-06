#include <RcppArmadillo.h>
#include "functions.h"

// [[Rcpp::export]]
arma::vec grad_(arma::vec& logPars, int G, 
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
  
  arma::vec likeas = arma::zeros<arma::vec>(logPars.n_elem);
  double likelam = 0.0;
  double likeb = 0.0;
  double likeq = 0.0;
  double liketau = 0.0;
  double likea2 = 0.0;
  double likeb2 = 0.0;
  double likec1 = 0.0;
  
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

  arma::vec gradient = arma::zeros<arma::vec>(logPars.n_elem);
  
  for(int i=0; i<m; i++){
    
    int mint_i = startSamplingPeriod[i];
    
    // int mint_i;
    // if(birthTimes[i]>1L){
    //   mint_i = birthTimes[i];
    // }else{
    //   mint_i = 1L;
    // }
    
    if(birthTimes[i] < startSamplingPeriod[i]){
      likea2 += DlogS_a2(ageMat(i, startSamplingPeriod[i] - birthTimes[i]), a2, b2);
      likeb2 += DlogS_b2(ageMat(i, startSamplingPeriod[i] - birthTimes[i]), a2, b2);
      likec1 += DlogS_c1(ageMat(i, startSamplingPeriod[i] - birthTimes[i]), c1);
    }
    
    for(int j=mint_i; j<lastObsAliveTimes[i]; j++){
      
      int g = SocGroup(i, j-1);
      
      // if(g==0L){
      //   Rcout << "g at grad = " << g << std::endl;
      //   Rcpp::stop("found g=0.");
      // }
      
      
      int z_t_1 = X(i,j-1);
      int z_t = X(i,j);
      
      if( ((z_t_1==0L) || (z_t_1==1L) || (z_t_1==3L)) && (z_t==9L) ){
        double age_ij = (double)(ageMat(i,j));
        double ptiOverqti = TrProbSurvive_(age_ij, a2, b2, c1, false) / 
          TrProbDeath_(age_ij, a2, b2, c1, false);
        likea2 -= ptiOverqti*Dlogpt_a2(age_ij, a2, b2);
        likeb2 -= ptiOverqti*Dlogpt_b2(age_ij, a2, b2);
        likec1 -= ptiOverqti*Dlogpt_c1(c1);
      }else if((z_t_1==0L) && (z_t==0L)){
        double mdenom = (double)(totalmPerGroup(g-1, j-1))/K;
        double inf_mgt = totalNumInfec(g-1, j-1)/(pow(mdenom, q));
        a = alpha_js[g-1];
        likeas[g-1] -= a;
        likelam -= a;
        likeb -= b*inf_mgt;
        likeq += b*inf_mgt*log(mdenom)*exp(ql)/(pow((1.0+exp(ql)), 2));
        double age_ij = (double)(ageMat(i,j));
        likea2 += Dlogpt_a2(age_ij, a2, b2);
        likeb2 += Dlogpt_b2(age_ij, a2, b2);
        likec1 += Dlogpt_c1(c1);
      }else if((z_t_1==1L) && (z_t==1L)){
        double age_ij = (double)(ageMat(i,j));
        likea2 += Dlogpt_a2(age_ij, a2, b2);
        likeb2 += Dlogpt_b2(age_ij, a2, b2);
        likec1 += Dlogpt_c1(c1);
      }else if((z_t_1==0L) && (z_t==3L)){
        double mdenom = (double)(totalmPerGroup(g-1, j-1))/K;
        double inf_mgt = totalNumInfec(g-1, j-1)/(pow(mdenom, q));
        a = alpha_js[g-1];
        double toBeExp = a+b*inf_mgt;
        if(toBeExp<1e-15){
          likeas[g-1] += 1.0;
          likelam += 1.0;
          if(totalNumInfec(g-1, j-1)==0L){
            likeb += 0.0;
          }else{
            likeb += 1.0;
          }
        }else{
          // double ratio = exp(-a-b*inf_mgt)/(1-exp(-a-b*inf_mgt));
          double ratio = exp( -a-b*inf_mgt - Rf_log1mexp(a+b*inf_mgt) );
          likeas[g-1] += a*ratio;
          likelam += a*ratio;
          likeb += b*ratio*inf_mgt;
          likeq -= b*ratio*inf_mgt*log(mdenom)*exp(ql)/(pow((1.0+exp(ql)), 2));
        }
        
        double age_ij = (double)(ageMat(i,j));
        likea2 += Dlogpt_a2(age_ij, a2, b2);
        likeb2 += Dlogpt_b2(age_ij, a2, b2);
        likec1 += Dlogpt_c1(c1);
      }else if((z_t_1==3L) && (z_t==3L)){
        liketau -= DerivErlangCDF(1, k, tau/((double)k))/(1-ErlangCDF(1, k, tau/((double)k)));
        double age_ij = (double)(ageMat(i,j));
        likea2 += Dlogpt_a2(age_ij, a2, b2);
        likeb2 += Dlogpt_b2(age_ij, a2, b2);
        likec1 += Dlogpt_c1(c1);
      }else if((z_t_1==3L) && (z_t==1L)){
        liketau += DerivErlangCDF(1, k, tau/((double)k))/ErlangCDF(1, k, tau/((double)k));
        double age_ij = (double)(ageMat(i,j));
        likea2 += Dlogpt_a2(age_ij, a2, b2);
        likeb2 += Dlogpt_b2(age_ij, a2, b2);
        likec1 += Dlogpt_c1(c1);
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
      likea2 += Dlogpt_a2(age_ij, a2, b2);
      likeb2 += Dlogpt_b2(age_ij, a2, b2);
      likec1 += Dlogpt_c1(c1);
    }
    
  }
  
  
  
  for(int g=0; g<G; g++){
    likeas[g] += 1.0 - exp(logPars[g]);
  }

  likelam += 1.0 - hp_lambda[1] * lambda;
  likeb += hp_beta[0] - hp_beta[1] * b;
  likeq += hp_q[0] - (hp_q[0]+hp_q[1])*q;
  liketau += hp_tau[0] - hp_tau[1] * tau;
  likea2 += hp_a2[0] - hp_a2[1] * a2;
  likeb2 += hp_b2[0] - hp_b2[1] * b2;
  likec1 += hp_c1[0] - hp_c1[1] * c1;


  for(int g=0; g<G; g++){
    gradient[g] = likeas[g];
  }
  gradient[G] = likelam;
  gradient[G+1] = likeb;
  gradient[G+2] = likeq;
  gradient[G+3] = liketau;
  gradient[G+4] = likea2;
  gradient[G+5] = likeb2;
  gradient[G+6] = likec1;
  
  return(gradient);
}
