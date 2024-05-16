#include <RcppArmadillo.h>
#include "funs_overcam.h"
#include "funs_cam.h"
#include "funs_fcam.h"
#include "funs_ficam.h"
#include "funs_san.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
Rcpp::List sample_overficam_burn(int nrep, // number of replications of the Gibbs sampler
                                 int burn,
                                 const arma::vec & y, // input data
                                 const arma::vec & group, // group assignment for each observation in the vector y
                                 int maxK, // maximum number of distributional clusters
                                 int maxL, // maximum number of observational clusters
                                 double m0, double tau0, // hyperparameters on the N-iG prior on the mean parameter, mu|sigma2 ~ N(m0, sigma2 / tau0)
                                 double lambda0, double gamma0, // hyperparameters on the N-iG prior on the variance parameter, 1/sigma2 ~ Gamma(lambda0, gamma0)
                                 bool fixed_alpha, bool fixed_beta, // do you want fixed alpha or beta?
                                 double alpha, double beta, // Dirichlet parameters if fixed
                                 double hyp_alpha1, double hyp_alpha2, // hyperparameters of alpha ( par of the distributional DP )
                                 double hyp_beta, // hyperparameter of the Gamma prior for the observational Dirichlet
                                 arma::vec mu_start, // starting point
                                 arma::vec sigma2_start,
                                 arma::vec M_start,
                                 arma::vec S_start,
                                 double alpha_start, double beta_start,
                                 double eps_beta, // MH step on beta
                                 bool progressbar
)
{
  int N = y.n_elem ;
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  
  
  // allocate output matrices
  arma::mat mu(maxL, nrep - burn, arma::fill::zeros) ; 
  arma::mat sigma2(maxL, nrep - burn, arma::fill::ones) ; 
  arma::mat out_S(G, nrep - burn) ; // distributional clusters
  arma::mat out_M(N, nrep - burn) ; // observational clusters
  arma::mat out_pi(maxK, nrep - burn, arma::fill::zeros) ;
  arma::cube out_omega(maxL, maxK, nrep - burn, arma::fill::zeros) ; 
  arma::vec out_alpha(nrep - burn) ; // DP parameter for the distributional weights 
  arma::vec out_beta(nrep - burn) ; // DP parameter for the observational weights 
  arma::vec out_maxK(nrep - burn) ; 

  // initialization
  arma::vec current_mu = mu_start ;
  arma::vec current_sigma2 = sigma2_start ;
  arma::vec current_M = M_start;
  arma::vec current_S = S_start;
  arma::vec current_pi(maxK) ; 
  current_pi.fill(1.0/(maxK));
  arma::mat current_omega(maxL,maxK) ; 
  current_omega.fill(1.0/(maxL));
  int current_maxK = maxK ;

  Rcpp::List out_params ;
  Rcpp::List tmp_obs_cl ;
  Rcpp::List weights_slice_sampler ;
  
  double current_beta = 0;
  double current_alpha = 0;
  
  if(fixed_alpha) { 
    out_alpha.fill(alpha) ; 
    current_alpha = alpha_start ;
  } else {
    current_alpha = alpha_start ;
  }
  
  if(fixed_beta) { 
    out_beta.fill(beta) ; 
    current_beta = beta_start;
  } else {
    current_beta = beta_start ;
  }
  
  arma::vec clusterD_long(N) ;
  for(int i = 0; i < N; i++) { clusterD_long(i) = current_S( group(i)) ; }
  arma::vec xi(maxK) ;
  for(int k = 0; k < maxK; k++) { xi(k) = fun_xi(0.5, k+1) ; }
  
  bool warningK = false ; 
  arma::vec iters_maxK(nrep) ;  // here changed
  int countK = 0 ;
  
  // progress bar
  bool display_progress = progressbar ;
  Progress p(nrep, display_progress) ;
  
  // START
  for(int iter = 0; iter < nrep ; iter++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    
    /*---------------------------------------------*/
      /*
      *  UPDATE MIXTURE WEIGHTS
    */
      
      
      weights_slice_sampler = ficam_weights_update_slice_sampler(group,
                                                                 current_S,
                                                                 current_alpha,
                                                                 xi,
                                                                 maxK) ;
      
      
      arma::vec tmp_pi = weights_slice_sampler["new_pi"] ;
      current_pi = tmp_pi ;
      
      int maxK_new = weights_slice_sampler["new_maxK"] ;
      current_maxK = maxK_new ;
      if(current_maxK >= maxK) { 
        warningK = true ;
        iters_maxK(countK) = (iter + 1) ;
        countK++ ;
        current_maxK = maxK-1 ; }    
      
    arma::vec u_D = weights_slice_sampler["u_D"] ;
    
    current_omega = dirichlet_sample_obs_weights(current_M,
                                                 clusterD_long,
                                                 current_beta*maxL,
                                                 current_maxK, maxL, // checked, added +1
                                                 maxK, maxL) ;
    
    /*---------------------------------------------*/
      
      
      /*---------------------------------------------*/
      /*
      *  UPDATE CLUSTER ASSIGNMENT
    */
      /* update distributional clusters S */
      current_S = slicedDP_sample_distr_cluster(group,
                                                        current_M,
                                                        current_pi, 
                                                        current_omega,
                                                        u_D, xi,
                                                        current_maxK) ;
    for(int i = 0; i < N; i++) { clusterD_long(i) = current_S( group(i)) ; }
    
    /* update observational clusters M */
      current_M = sample_obs_cluster(y, 
                                     clusterD_long,
                                     current_omega,
                                     current_maxK, maxL,
                                     current_mu, 
                                     current_sigma2) ;
    
    /*---------------------------------------------*/
      
      
      /*---------------------------------------------*/
      /*
      *  UPDATE PARAMETERS
    */
      out_params = sample_model_parameters(y, current_M,
                                           maxL,
                                           m0, tau0, lambda0, gamma0) ;
    
    arma::vec tmp_mu = out_params["out_mu"] ;
    current_mu = tmp_mu ; 
    
    arma::vec tmp_sigma2 = out_params["out_sigma2"] ;
    current_sigma2 = tmp_sigma2 ; 
    
    /*---------------------------------------------*/
      
      /*---------------------------------------------*/
      /*
      *  UPDATE DIRICHLET HYPER-PARAMETERS
    */
      /* sample alpha */
      if(!fixed_alpha) {
        current_alpha = sample_alpha(current_alpha,
                                     hyp_alpha1, hyp_alpha2,
                                     G, current_S) ;
      }
      
    /* sample beta */
      if(!fixed_beta) {
        current_beta = overcam_MH_beta(current_beta, eps_beta, 
                                       current_omega, hyp_beta) ;
      }
      

    /*---------------------------------------------*/
      
      // saving
      
      if(iter >= burn){
        int ind = iter-burn;
        //Rcpp::Rcout << ind << "\n";
        
        mu.col(ind) = current_mu;
        sigma2.col(ind) = current_sigma2;
        //Rcpp::Rcout << "a";
        
        out_alpha(ind) = current_alpha;  
        out_beta(ind)  = current_beta;  
        //Rcpp::Rcout << "c";
        
        out_pi.col(ind) = current_pi;
        out_omega.slice(ind) = current_omega;
        //      Rcpp::Rcout << "d";
        
        out_S.col(ind) = current_S;
        out_M.col(ind) = current_M;
        //      Rcpp::Rcout << "e";
        
      }
      
      // END
  }
  
  Rcpp::List warnings ;
  if(warningK) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxK") = iters_maxK.head(countK) );
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu.t(),
                            Rcpp::Named("sigma2") = sigma2.t(),
                            Rcpp::Named("obs_cluster") = out_M.t(),
                            Rcpp::Named("distr_cluster") = out_S.t(),
                            Rcpp::Named("pi") = out_pi.t(),
                            Rcpp::Named("omega") = out_omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("maxK") = out_maxK,
                            Rcpp::Named("warnings") = warnings);
}
