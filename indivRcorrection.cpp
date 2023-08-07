#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::cube CalcIndivRcorretion(arma::mat NAmatrix,
                         arma::imat SocGroup, 
                         arma::imat infTimes,
                         arma::imat infectivityTimes,
                         arma::imat deathTimes,
                         arma::icube nInfByGroup,
                         arma::icube nTotByGroup,
                         arma::mat alphas, 
                         arma::vec beta,
                         arma::vec q_,
                         double K){
  
  int N = alphas.n_rows;
  int G = alphas.n_cols;
  int m = SocGroup.n_rows;
  int maxt = SocGroup.n_cols;
  
  arma::cube relRatesByIndivIter = arma::zeros<arma::cube>(m, maxt, N);
  
  arma::mat totalRates = arma::zeros<arma::mat>(G, maxt);
  
  arma::mat relRatesByIndiv = arma::zeros<arma::mat>(m, maxt);
  
  
  for(int iter=0; iter<N; iter++){
    
    arma::vec as = (alphas.row(iter)).t();
    double b = beta[iter];
    double q = q_[iter];
    double r_gt;
    double I;
    double pop;
    
    totalRates.fill(0.0);
    
    for(int g=0; g<G; g++){
      for(int t=0; t<maxt; t++){
        I = (double)(nInfByGroup(g,t,iter));
        if(nTotByGroup(g,t,iter)>0L){
          pop = (double)(nTotByGroup(g,t,iter));
          r_gt = as[g] + b*I/(pow(pop/K, q));
          // r_gt = exp(-b/pow(pop/K, q))*(1-exp(-as[g])) + 
          //   I*(1 - exp(-b/pow(pop/K, q)))*exp(-as[g]);
          totalRates(g,t) = r_gt;
        }
      }
    }
    
    // relRatesByIndiv.fill(-10.0);
    
    relRatesByIndiv = NAmatrix;
    
    double R_it;
    int g;
    double contrib_it;
    int tmin;    // -1 is used because of c++ syntax
    int tmax;    // -1 is used because of (non-included) upper bound (death time -1)
    for(int i=0; i<m; i++){
      
      // badger must be infected at some point
      // badger must have died at some point
      if((infTimes(i,iter)!=-10L) && (deathTimes(i,iter)!=-10L)){
        
        
        if(infectivityTimes(i,iter)==-10L){  // if never infectious
          
          tmin = infTimes(i,iter) - 1L;
          tmax = deathTimes(i,iter) - 1L;
          
          for(int t=tmin; t<tmax; t++){
            relRatesByIndiv(i,t) = 0.0;
          }
          
        }else{   // if infectious at some time
          
          tmin = infTimes(i,iter) - 1L;
          tmax = infectivityTimes(i,iter) - 1L;
          
          for(int t=tmin; t<tmax; t++){
            relRatesByIndiv(i,t) = 0.0;
          }
          
          tmin = infectivityTimes(i,iter) - 1L;
          tmax = deathTimes(i,iter) - 1L;
          
          for(int t=tmin; t<tmax; t++){
            
            g = SocGroup(i,t) - 1L;
            
            arma::uvec whichSEeventInGroupg = arma::find(
              (infTimes.col(iter) == t+1) && (SocGroup.col(t)==g+1L));
            
            int numEvents = whichSEeventInGroupg.n_elem;
            
            if(numEvents>0){
              pop = (double)(nTotByGroup(g,t,iter));
              if(nTotByGroup(g,t,iter)==0L){
                Rcpp::stop("pop = 0");
              }
              
              
              contrib_it = b/(pow(pop/K, q));
              // contrib_it = (1 - exp(-b/pow(pop/K, q)))*exp(-as[g]);
              
              R_it = (double)(numEvents) * contrib_it / totalRates(g,t);
              relRatesByIndiv(i,t) = R_it;
            }else{
              relRatesByIndiv(i,t) = 0.0;
            }
            
          } // end loop over infectious period
          
        }
        
        
      } // end if (infected and dead)
      
    } // end loop over individuals
    
    
    relRatesByIndivIter.slice(iter) = relRatesByIndiv;
    
  } // end loop over iterations
  
  return(relRatesByIndivIter);
  
  // List outList(1);
  // outList[0] = relRatesByIndivIter;
  // return(outList);
  
}
