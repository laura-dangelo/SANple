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
Rcpp::List sample_overficam_arma(int nrep, // number of replications of the Gibbs sampler
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
  arma::mat mu(maxL, nrep, arma::fill::zeros) ;
  arma::mat sigma2(maxL, nrep, arma::fill::ones) ; 
  arma::mat out_S(G, nrep) ;
  arma::mat out_M(N, nrep) ;
  arma::mat pi(maxK, nrep, arma::fill::zeros) ;
  arma::cube omega(maxL, maxK, nrep, arma::fill::zeros) ;
  arma::vec out_alpha(nrep) ; // DP parameter for the distributional weights
  arma::vec out_beta(nrep) ; // Dirichlet parameter for the observational weights
  arma::vec out_maxK(nrep) ;

  // initialization
  mu.col(0) = mu_start ;
  sigma2.col(0) = sigma2_start ;
  out_M.col(0) = M_start ;
  out_S.col(0) = S_start ;
  out_maxK(0) = maxK ;
  pi.col(0).fill(1/(maxK*1.0)) ;
  omega.slice(0).fill(1/(maxL*1.0)) ;

  Rcpp::List out_params ;
  Rcpp::List tmp_obs_cl ;
  Rcpp::List weights_slice_sampler ;

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
  arma::vec xi(maxK) ;
  for(int k = 0; k < maxK; k++) { xi(k) = fun_xi(0.5, k+1) ; }
  
  bool warningK = false ; arma::vec iters_maxK(nrep-1) ; int countK = 0 ;
  
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

    weights_slice_sampler = ficam_weights_update_slice_sampler(group,
                                                               out_S.col(iter),
                                                               out_alpha(iter),
                                                               xi,
                                                               maxK) ;

    arma::vec tmp_pi = weights_slice_sampler["new_pi"] ;
    pi.col(iter+1) = tmp_pi ;

    int maxK_new = weights_slice_sampler["new_maxK"] ;
    out_maxK(iter+1) = maxK_new ;
    if(out_maxK(iter + 1) >= maxK) { 
      warningK = true ;
      iters_maxK(countK) = (iter + 1) ;
      countK++ ;
      out_maxK(iter + 1) = maxK-1 ; }
    
    arma::vec u_D = weights_slice_sampler["u_D"] ;

    omega.slice(iter+1) = dirichlet_sample_obs_weights(out_M.col(iter),
                                                       clusterD_long,
                                                       out_beta(iter)*maxL,
                                                       out_maxK(iter+1), maxL,
                                                       maxK, maxL) ;
    
    /*---------------------------------------------*/
    

    /*---------------------------------------------*/
    /*
     *  UPDATE CLUSTER ASSIGNMENT
     */
    /* update distributional clusters S */
    out_S.col(iter+1) = slicedDP_sample_distr_cluster(y, group,
                                                      out_M.col(iter),
                                                      pi.col(iter+1), omega.slice(iter+1),
                                                      u_D, xi,
                                                      out_maxK(iter+1)) ;
    for(int i = 0; i < N; i++) { clusterD_long(i) = out_S( group(i), iter+1 ) ; }
    
    /* update observational clusters M */
    out_M.col(iter+1) = sample_obs_cluster(y, 
                                           clusterD_long,
                                           omega.slice(iter+1),
                                           out_maxK(iter+1), maxL,
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
      out_alpha(iter+1) = sample_alpha(out_alpha(iter),
                                      hyp_alpha1, hyp_alpha2,
                                      G, out_S.col(iter+1)) ;
    }
    
    /* sample beta */
    if(!fixed_beta) {
      out_beta(iter+1) = overcam_MH_beta(out_beta(iter), eps_beta, omega.slice(iter+1), hyp_beta) ;
    }
    
    /*---------------------------------------------*/
    

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
                            Rcpp::Named("pi") = pi.t(),
                            Rcpp::Named("omega") = omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("maxK") = out_maxK,
                            Rcpp::Named("warnings") = warnings);
}
