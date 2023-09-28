#include "funs_ficam.h"
// compute cluster probabilities and sample new cluster allocation
// the function outputs:
// distributional cluster allocation variables (S)
// distributional prob pi


Rcpp::List ficam_weights_update_slice_sampler(const arma::vec& group, 
                                        arma::vec S_iter, 
                                        double alpha,
                                        arma::vec xi,
                                        int maxK)
{
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  arma::vec out_pi(maxK, arma::fill::zeros) ;
  
  int maxK_new ; 
  arma::vec u_D(G) ;
  
  int a_k ; double b_k ;
  
  /*---------------------------------------------*/
  /*
  *  SAMPLE SLICE - uniform r.v.
  */
  /* distributional */
  for(int j = 0; j < G; j++)
  {
    u_D(j) = R::runif( 0, xi(S_iter(j)) ) ;
  }
  int maxK_slice_tmp = compute_trunc(u_D.min(), 0.5) + 1 ;
  maxK_new = std::min(maxK, maxK_slice_tmp) ;
  
  /*---------------------------------------------*/
    
    
  /*---------------------------------------------*/
  /*
  *  UPDATE MIXTURE WEIGHTS
  */
  /* sample distributional probabilities pi */
  /* (pi_1,...,pi_maxK) | . ~ DP */
    
  arma::vec v_k(maxK_new) ;
  arma::vec pi_k(maxK_new) ; 
  for(int k = 0; k < maxK_new ; k++)
  {
    arma::uvec ind_k = find(S_iter == k) ;
    arma::uvec ind_mk = find(S_iter > k) ;

    a_k = 1 + ind_k.n_elem ;
    b_k = alpha + ind_mk.n_elem ;
    v_k(k) = R::rbeta(a_k, b_k) ;
  }
  pi_k = stick_breaking( v_k ) ;
  out_pi(arma::span(0, maxK_new-1)) = pi_k ;
  
  
  /*---------------------------------------------*/
    
  return Rcpp::List::create(Rcpp::Named("new_maxK") = maxK_new,
                            Rcpp::Named("new_pi") = out_pi,
                            Rcpp::Named("u_D") = u_D);
  
}
