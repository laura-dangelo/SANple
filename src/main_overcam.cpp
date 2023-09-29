#include <RcppArmadillo.h>
#include "funs_overcam.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
Rcpp::List sample_overcam_arma(int nrep, // number of replications of the Gibbs sampler
                               const arma::vec & y, // input data
                               const arma::vec & group, // group assignment for each observation in the vector y
                               int maxK, // maximum number of distributional clusters
                               int maxL, // maximum number of observational clusters  
                               double m0, double tau0, // hyperparameters on the N-iG prior on the mean parameter, mu|sigma2 ~ N(m0, sigma2 / tau0)
                               double lambda0, double gamma0, // hyperparameters on the N-iG prior on the variance parameter, 1/sigma2 ~ Gamma(lambda0, gamma0)
                               bool fixed_alpha, bool fixed_beta, // do you want fixed alpha or beta?
                               double alpha, double beta, // Dirichlet parameters if fixed
                               double hyp_alpha, // hyperparameter of the Gamma prior for the distributional Dirichlet
                               double hyp_beta, // hyperparameter of the Gamma prior for the observational Dirichlet
                               arma::vec mu_start, // starting point 
                               arma::vec sigma2_start,
                               arma::vec M_start,
                               arma::vec S_start,
                               double alpha_start,
                               double beta_start,
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
  
  // initialization
  mu.col(0) = mu_start ;
  sigma2.col(0) = sigma2_start ;
  out_M.col(0) = M_start ;
  out_S.col(0) = S_start ;
  
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
  
  // other auxiliary quantities
  Rcpp::List out_params ;
  arma::vec clusterD_long(N) ; 
  for(int i = 0; i < N; i++) { clusterD_long(i) = out_S( group(i), 0 ) ; }
  
  arma::vec iters_maxL(nrep-1);
  arma::vec iters_maxK(nrep-1);
  
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
     *  UPDATE MIXTURE WEIGHTS
     */
    /* sample distributional probabilities pi */
    pi.col(iter+1) = dirichlet_sample_distr_weights(out_S.col(iter), 
                                                    out_alpha(iter)*maxK, 
                                                    maxK, 
                                                    maxK) ;
    
    /* for each distributional cluster k, sample observational probabilities omega */
    omega.slice(iter+1) = dirichlet_sample_obs_weights(out_M.col(iter), 
                                                       clusterD_long,
                                                       out_beta(iter)*maxL,
                                                       maxK, maxL,
                                                       maxK, maxL) ;
    /*---------------------------------------------*/
    
    
    /*---------------------------------------------*/
    /*
     *  UPDATE CLUSTER ASSIGNMENT
     */
    
    /* update distributional clusters S */
    out_S.col(iter+1) = sample_distr_cluster(y, group,
                                             pi.col(iter+1), omega.slice(iter+1),
                                             maxK, maxL,
                                             mu.col(iter), sigma2.col(iter)) ;
    
    for(int i = 0; i < N; i++) { clusterD_long(i) = out_S( group(i), iter+1 ) ; }
    
    /* update observational clusters M */
    out_M.col(iter+1) = sample_obs_cluster(y, 
                                           clusterD_long,
                                           omega.slice(iter+1),
                                           maxK, maxL,
                                           mu.col(iter), sigma2.col(iter)) ;
    
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
     *  UPDATE DIRICHLET HYPER-PARAMETERS 
     */
    /* sample alpha */
    if(!fixed_alpha) {
      out_alpha(iter+1) = overcam_MH_alpha(out_alpha(iter), eps_alpha, pi.col(iter+1), hyp_alpha) ;
    }
    
    /* sample beta */
    if(!fixed_beta) {
      out_beta(iter+1) = overcam_MH_beta(out_beta(iter), eps_beta, omega.slice(iter+1), hyp_beta) ;
    }
    /*---------------------------------------------*/
    
    // END
  }
  
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu.t(),
                            Rcpp::Named("sigma2") = sigma2.t(),
                            Rcpp::Named("obs_cluster") = out_M.t(),
                            Rcpp::Named("distr_cluster") = out_S.t(),
                            Rcpp::Named("pi") = pi.t(),
                            Rcpp::Named("omega") = omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta);
}
