#include <RcppArmadillo.h>
#include "funs_fcam.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


// [[Rcpp::export]]
Rcpp::List sample_fcam_burn(int nrep, // number of replications of the Gibbs sampler
                            int burn,
                            const arma::vec & y, // input data
                            const arma::vec & group, // group assignment for each observation in the vector y
                            int maxK, // maximum number of distributional clusters
                            int maxL, // maximum number of observational clusters  
                            double m0, double tau0, // hyperparameters on the N-iG prior on the mean parameter, mu|sigma2 ~ N(m0, sigma2 / tau0)
                            double lambda0, double gamma0, // hyperparameters on the N-iG prior on the variance parameter, 1/sigma2 ~ Gamma(lambda0, gamma0)
                            bool fixed_alpha, bool fixed_beta, // do you want fixed alpha or beta?
                              double alpha, double beta, // Dirichlet parameters if fixed
                            double hyp_alpha1, double hyp_alpha2,// hyperparameter of the gamma prior for the distributional Dirichlet
                            double hyp_beta1, double hyp_beta2, // hyperparameter of the gamma prior for the observational Dirichlet
                            arma::vec mu_start, // starting point 
                            arma::vec sigma2_start,
                            arma::vec M_start,
                            arma::vec S_start,
                            double alpha_start, double beta_start,
                            double eps_alpha, // MH step on alpha
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
  arma::vec out_K(nrep - burn) ; 
  arma::vec out_L(nrep - burn) ;
  arma::vec out_Kplus(nrep-burn) ;
  arma::vec out_Lplus(nrep-burn) ;
  
  // initialization
  arma::vec current_mu = mu_start ;
  arma::vec current_sigma2 = sigma2_start ;
  arma::vec current_M = M_start ;
  arma::vec current_S = S_start ;
  int current_K = maxK-1 ;
  int current_L = maxL-1 ;
  int current_Kplus = current_S.max() + 1 ;
  int current_Lplus = current_M.max() + 1 ;
  arma::vec current_pi = arma::ones(maxK)/maxK ;
  arma::mat current_omega = arma::ones(maxL, maxK) / maxL ;
  
  Rcpp::List out_params ;
  Rcpp::List tmp_obs_cl ;
  bool warningL = false ; arma::vec iters_maxL(nrep) ; int countL = 0 ;
  bool warningK = false ; arma::vec iters_maxK(nrep) ; int countK = 0 ;
  
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
  
  
  
  // auxiliary quantities
  arma::vec clusterD_long(N) ; 
  for(int i = 0; i < N; i++) { clusterD_long(i) = S_start( group(i) ) ; }

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
      *  UPDATE CLUSTER ASSIGNMENT
    */
      /* update distributional clusters S */
      current_S = fcam_sample_distr_cluster(y, group,
                                            current_pi, 
                                            current_omega,
                                            current_M,
                                            current_K) ;
    
    current_Kplus = current_S.max() + 1 ;
    for(int i = 0; i < N; i++) { clusterD_long(i) = current_S( group(i) ) ; }
    
    /* update observational clusters M */
      current_M = fcam_sample_obs_cluster(y, 
                                                  clusterD_long,
                                                  current_omega,
                                                  current_K, 
                                                  current_L,
                                                  current_mu, 
                                                  current_sigma2) ;
    
    // out_M.col(iter+1) = out_M.col(iter) ;
    current_Lplus = current_M.max() + 1;
    
    
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
      *  UPDATE THE NUMBER OF COMPONENTS
    */
      
      /* update distributional components out_K */
      current_K = fcam_sample_K(maxK,
                                      current_Kplus,
                                      1, 4, 3,
                                      current_alpha,
                                      current_S) ;
    if(current_K >= maxK) { 
      warningK = true ;
      iters_maxK(countK) = (iter + 1) ;
      countK++ ;
      current_K = maxK-1 ; }
    
    /* update observational clusters out_L */
      current_L = fcam_sample_K(maxL,
                                      current_Lplus,
                                      1, 4, 3,
                                      current_beta,
                                      current_M) ;
    if(current_L >= maxL) { 
      warningL = true ;
      iters_maxL(countL) = (iter + 1) ;
      countL++ ;
      current_L = maxL-1 ; }
    
    current_mu(arma::span(current_L, maxL-1)) = arma::zeros(maxL - current_L) ;
    current_sigma2(arma::span(current_L, maxL-1)) = arma::zeros(maxL -current_L) ;
    
    /*---------------------------------------------*/
      
      /*---------------------------------------------*/
      /*
      *  UPDATE DIRICHLET HYPER-PARAMETERS 
    */
      /* sample alpha */
      if(!fixed_alpha) {
        current_alpha = fcam_MH_alpha(current_alpha, eps_alpha,
                                          hyp_alpha1, hyp_alpha2,
                                          current_Kplus, current_K,
                                          current_S) ;
      }
    
    /* sample beta */
      if(!fixed_beta) {
        current_beta = fcam_MH_alpha(current_beta, eps_beta,
                                     hyp_beta1, hyp_beta2,
                                     current_Lplus, current_L,
                                     current_M) ;
      }
    /*---------------------------------------------*/
      
      
      /*---------------------------------------------*/
      /*
      *  UPDATE MIXTURE WEIGHTS
    */
      current_pi = dirichlet_sample_distr_weights(current_S,
                                                  current_alpha,
                                                  current_K, maxK) ;
    
    current_omega = dirichlet_sample_obs_weights(current_M,
                                                 clusterD_long,
                                                 current_beta,
                                                 current_K+1, current_L,
                                                 maxK, maxL) ;
    
    
    // pi.col(iter+1) = pi.col(iter) ;
    // omega.slice(iter+1) = omega.slice(iter) ;
    /*---------------------------------------------*/
      
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
       
       out_Kplus(ind) = current_Lplus;
       out_Lplus(ind) = current_Kplus;
       out_L(ind) = current_L;
       out_K(ind) = current_K;
      }
      
      
      
      // END
  }
  
  Rcpp::List warnings ;
  if(warningK & !(warningL)) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxK") = iters_maxK.head(countK) );
  }
  if(warningL & !(warningK)) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxL") = iters_maxL.head(countL) );
  }
  if(warningL & warningK) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxK") = iters_maxK.head(countK),
                                  Rcpp::Named("top_maxL") = iters_maxL.head(countL) );
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu.t(),
                            Rcpp::Named("sigma2") = sigma2.t(),
                            Rcpp::Named("obs_cluster") = out_M.t(),
                            Rcpp::Named("distr_cluster") = out_S.t(),
                            Rcpp::Named("pi") = out_pi.t(),
                            Rcpp::Named("omega") = out_omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("K_plus") = out_Kplus,
                            Rcpp::Named("L_plus") = out_Lplus,
                            Rcpp::Named("K") = out_K,
                            Rcpp::Named("L") = out_L,
                            Rcpp::Named("warnings") = warnings);
}
