#include <RcppArmadillo.h>
#include "funs_fcam.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

 
// [[Rcpp::export]]
Rcpp::List sample_fcam_arma(int nrep, // number of replications of the Gibbs sampler
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
  arma::mat mu(maxL, nrep, arma::fill::zeros) ; 
  arma::mat sigma2(maxL, nrep, arma::fill::ones) ; 
  arma::mat out_S(G, nrep) ;
  arma::mat out_M(N, nrep) ;
  arma::mat pi(maxK, nrep, arma::fill::zeros) ;
  arma::cube omega(maxL, maxK, nrep, arma::fill::zeros) ; 
  arma::vec out_alpha(nrep) ; // Dirichlet parameter for the distributional weights 
  arma::vec out_beta(nrep) ; // Dirichlet parameter for the observational weights 
  arma::vec out_K(nrep) ;
  arma::vec out_L(nrep) ;
  arma::vec out_Kplus(nrep) ;
  arma::vec out_Lplus(nrep) ;
  
  // initialization
  mu.col(0) = mu_start ;
  sigma2.col(0) = sigma2_start ;
  out_M.col(0) = M_start ;
  out_S.col(0) = S_start ;
  out_K(0) = maxK-1 ;
  out_L(0) = maxL-1 ;
  out_Kplus(0) = out_S.col(0).max() + 1 ;
  out_Lplus(0) = out_M.col(0).max() + 1 ;
  pi.col(0) = arma::ones(maxK)/maxK ;
  omega.slice(0) = arma::ones(maxL, maxK) / maxL ;
  
  Rcpp::List out_params ;
  Rcpp::List tmp_obs_cl ;
  bool warningL = false ; arma::vec iters_maxL(nrep-1) ; int countL = 0 ;
  bool warningK = false ; arma::vec iters_maxK(nrep-1) ; int countK = 0 ;
  
  if(fixed_alpha) { 
    out_alpha.fill(alpha) ; 
  } else {
    out_alpha(0) = alpha_start ;
  }
  
  if(fixed_beta) { 
    out_beta.fill(beta) ; 
  } else {
    out_beta(0) = beta_start ; 
  }
  
  arma::vec clusterD_long(N) ; 
  for(int i = 0; i < N; i++) { clusterD_long(i) = out_S( group(i), 0 ) ; }
  
  // progress bar 
  bool display_progress = progressbar ;
  Progress p(nrep, display_progress) ;
  
  // START
  for(int iter = 0; iter < nrep -1 ; iter++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    
    /*---------------------------------------------*/
    /*
     *  UPDATE CLUSTER ASSIGNMENT
     */
    /* update distributional clusters S */
    out_S.col(iter+1) = fcam_sample_distr_cluster(y, group,
              pi.col(iter), omega.slice(iter),
              out_M.col(iter),
              out_K(iter)) ;
    
    out_Kplus(iter+1) = out_S.col(iter+1).max() + 1 ;
    for(int i = 0; i < N; i++) { clusterD_long(i) = out_S( group(i), iter+1 ) ; }
    
    /* update observational clusters M */
    out_M.col(iter+1) = fcam_sample_obs_cluster(y, 
              clusterD_long,
              omega.slice(iter),
              out_K(iter), out_L(iter),
              mu.col(iter), sigma2.col(iter)) ;
    
    // out_M.col(iter+1) = out_M.col(iter) ;
    out_Lplus(iter+1) = out_M.col(iter+1).max() + 1;
    
    
    /*---------------------------------------------*/
    
    
    /*---------------------------------------------*/
    /*
     *  UPDATE PARAMETERS
     */
    out_params = sample_model_parameters(y, out_M.col(iter+1),
                                         maxL,
                                         m0, tau0, lambda0, gamma0) ;
    
    arma::vec tmp_mu = out_params["out_mu"] ;
    mu.col(iter+1) = tmp_mu ; 
    
    arma::vec tmp_sigma2 = out_params["out_sigma2"] ;
    sigma2.col(iter+1) = tmp_sigma2 ; 
    
    /*---------------------------------------------*/
    
    
    /*---------------------------------------------*/
    /*
     *  UPDATE THE NUMBER OF COMPONENTS
     */
    
    /* update distributional components out_K */
    out_K(iter + 1) = fcam_sample_K(maxK,
                                    out_Kplus(iter+1),
                                    1, 4, 3,
                                    out_alpha(iter),
                                    out_S.col(iter + 1)) ;
    if(out_K(iter + 1) >= maxK) { 
      warningK = true ;
      iters_maxK(countK) = (iter + 1) ;
      countK++ ;
      out_K(iter + 1) = maxK-1 ; }
    
    /* update observational clusters out_L */
    out_L(iter + 1) = fcam_sample_K(maxL,
                                     out_Lplus(iter+1),
                                     1, 4, 3,
                                     out_beta(iter),
                                     out_M.col(iter + 1)) ;
    if(out_L(iter + 1) >= maxL) { 
      warningL = true ;
      iters_maxL(countL) = (iter + 1) ;
      countL++ ;
      out_L(iter + 1) = maxL-1 ; }
    
    mu(arma::span(out_L(iter+1), maxL-1), arma::span(iter+1)) = arma::zeros(maxL - out_L(iter+1)) ;
    sigma2(arma::span(out_L(iter+1), maxL-1), iter+1) = arma::zeros(maxL - out_L(iter+1)) ;

    /*---------------------------------------------*/
    
    /*---------------------------------------------*/
    /*
     *  UPDATE DIRICHLET HYPER-PARAMETERS 
     */
    /* sample alpha */
    if(!fixed_alpha) {
      out_alpha(iter+1) = fcam_MH_alpha(out_alpha(iter), eps_alpha,
                hyp_alpha1, hyp_alpha2,
                out_Kplus(iter+1), out_K(iter+1),
                out_S.col(iter+1)) ;
    }
    
    /* sample beta */
    if(!fixed_beta) {
      out_beta(iter+1) = fcam_MH_alpha(out_beta(iter), eps_beta,
               hyp_beta1, hyp_beta2,
               out_Lplus(iter+1), out_L(iter+1),
               out_M.col(iter+1)) ;
    }
    /*---------------------------------------------*/
    
    
    /*---------------------------------------------*/
    /*
     *  UPDATE MIXTURE WEIGHTS
     */
    pi.col(iter+1) = dirichlet_sample_distr_weights(out_S.col(iter+1),
                                                    out_alpha(iter+1),
                                                    out_K(iter+1), maxK) ;

    omega.slice(iter+1) = dirichlet_sample_obs_weights(out_M.col(iter+1),
                                                       clusterD_long,
                                                       out_beta(iter+1),
                                                       out_K(iter+1)+1, out_L(iter+1),
                                                       maxK, maxL) ;

            
    // pi.col(iter+1) = pi.col(iter) ;
    // omega.slice(iter+1) = omega.slice(iter) ;
    /*---------------------------------------------*/
    
  
     
   
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
                            Rcpp::Named("pi") = pi.t(),
                            Rcpp::Named("omega") = omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("K_plus") = out_Kplus,
                            Rcpp::Named("L_plus") = out_Lplus,
                            Rcpp::Named("K") = out_K,
                            Rcpp::Named("L") = out_L,
                            Rcpp::Named("warnings") = warnings);
}
