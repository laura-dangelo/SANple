#include "funs_fcam.h"


arma::vec fcam_sample_distr_cluster(const arma::vec& y, const arma::vec& group,
                                    arma::vec pi, arma::mat omega,
                                    arma::vec M_iter,
                                    int K_iter)
{
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  
  arma::vec distr_cluster_id = arma::linspace(0, K_iter-1, K_iter) ;
  arma::vec probD(K_iter) ;
  
  arma::vec out(G) ;
  // S_j is categorical, with Pr( S_j = k | - ) = pi_k * ( om_1k^#(M_ij = 1) ... om_Lk^#(M_ij = L) ) 
  // or equivalently,  = pi_k * (om_{M_1j,k} * ... * om_{M_Nj,k})
  // hence the log is: log(pi_k) = log(pi_k) + sum_{i=1}^{Nj} log( om_{M_ij,k} ) 
  for(int j = 0; j < G; j++)
  {
    double sumdens ;
    arma::uvec ind_group_j = find(group == j) ;  // restrict to observations in the j-th population
    arma::vec mixdens(ind_group_j.n_elem) ; // empty vector of length N_j
    
    for(int k = 0; k < K_iter; k++)  // I have to compute the prob for each k = 1,..., maxK
    {
      mixdens.fill(0) ;
      for(int i = 0; i < ind_group_j.n_elem; i++) { mixdens(i) = log( omega(M_iter(ind_group_j(i)), k) ) ; }
      sumdens = arma::accu(mixdens) ;
      if(!arma::is_finite(sumdens)) { sumdens = log(0) ;}
      probD(k) =  log( pi(k) ) + sumdens ;
    }
    double tmp = max(probD) ;
    for(int k = 0; k < K_iter; k++) {  probD(k) = exp( probD(k) - tmp ) ; }
    out(j) = sample_i(distr_cluster_id, probD) ;
  }
  
  // relabel the clusters so that the first Kplus labels are used 
  out = relabel_arma(out) ;
  
  return(out) ;
}




arma::vec fcam_sample_obs_cluster(const arma::vec& y, 
                                  arma::vec clusterD_long,
                                  arma::mat omega,
                                  int K_iter, int L_iter,
                                  arma::vec mu, arma::vec sigma2)
{
  arma::vec obs_cluster_id = arma::linspace(0, L_iter-1, L_iter) ;
  int N = y.n_elem ;
  arma::vec probO(L_iter) ;
  
  arma::vec out(N) ;
  
  // M_ij is categorical, with Pr(M_ij = l | S_j = k) = om_lk * Pr(y_ij | mu_l)
  for(int i = 0; i < N; i++) 
  {
    for(int l = 0; l < L_iter; l++) {
      probO(l) = log( omega( l, clusterD_long(i) ) ) + R::dnorm(y(i), mu(l), std::sqrt(sigma2(l)), true) ;
    }
    probO =  probO - max(probO) ;
    probO = exp(probO) ;
    if(arma::accu(probO) == 0)  { Rcpp::Rcout << "all 0 in sampling M" ; probO = arma::ones(L_iter)/L_iter ;}
    out(i) = sample_i(obs_cluster_id, probO) ;
  }
  
  // relabel the clusters so that the first Lplus labels are used
  out = relabel_arma(out) ;
  
  return(out);
}





double lbeta(double a, double b)
{
  double out = 0;
  out = lgamma(a) + lgamma(b) - lgamma(a+b) ;
  return(out) ;
}

double logprior_maxx(int maxx, double hyp_max1,
                     double hyp_max2, double hyp_max3)
{
  double out ;
  out = lgamma( hyp_max1 + maxx - 1 ) + lbeta( hyp_max1 + hyp_max2, maxx - 1 + hyp_max3) -
    lgamma( hyp_max1 ) - lgamma( maxx ) - lbeta( hyp_max2, hyp_max3 ) ;
  return(out) ;
}

// log-posterior on max components 
double logpost_maxx(int Kiter, 
                    double hyp_max1, 
                    double hyp_max2, 
                    double hyp_max3, 
                    double alpha,
                    arma::vec& cluster,
                    int Kplus)
{
  double out ;
  
  double tmp = 0 ;
  for(int h = 0; h < Kplus - 1 ; h ++) 
  {
    arma::uvec ind = find( cluster == h ) ;
    tmp = tmp + lgamma( ind.n_elem + alpha/Kiter ) - lgamma( 1 + alpha/Kiter ) ;
  }
  
  out = logprior_maxx(Kiter, hyp_max1, hyp_max2, hyp_max3) + 
        Kplus * log(alpha) + lgamma(Kiter + 1) - Kplus * log(Kiter) - 
        lgamma(Kiter - Kplus + 1) + tmp ; 
  
  return(out) ;
}

int fcam_sample_K(int maxK,
                  int Kplus,
                  double hyp_max1, 
                  double hyp_max2, 
                  double hyp_max3,
                  double d_par,
                  arma::vec cluster)
{
  int newK ;
  arma::vec K_id = arma::linspace( Kplus, maxK, maxK - Kplus + 1 ) ;
  
  arma::vec prob2(K_id.n_elem) ;
  arma::vec prob(K_id.n_elem) ;
  for(int k = 0; k < K_id.n_elem ; k++) {
    prob(k) = logpost_maxx( Kplus + k , 
                            hyp_max1, hyp_max2, hyp_max3, 
                            d_par, cluster, Kplus ) ;
  }
  double minprob = prob.min() ;
  for(int k = 0; k < K_id.n_elem ; k++)  {
    prob2(k) = exp(prob(k) - minprob) ;
  }
  
  newK = sample_i(K_id, prob2) ;
  
  return(newK) ;
}



// prior on alpha and beta
double fcam_logprior_alpha(double alpha, 
                     double hyp_alpha1, double hyp_alpha2)
{
  double out ;
  // out = R::df(alpha, hyp_alpha1, hyp_alpha2, true) ;
  out = R::dgamma(alpha, hyp_alpha1, 1/hyp_alpha2, true) ;
  return(out) ;
}

// log-posterior for MH step on alpha and beta
double fcam_logpost_alpha(double alpha, 
                          double hyp_alpha1, double hyp_alpha2,
                          arma::vec cluster,
                          int Kplus,
                          int Kiter)
{
  double out ;
  int N = cluster.n_elem ;
  
  double tmp = 0 ;
  for(int h = 0; h < Kplus - 1 ; h ++) 
  {
    arma::uvec ind = find( cluster == h ) ;
    tmp = tmp + lgamma( ind.n_elem + alpha/Kiter ) - lgamma( 1 + alpha/Kiter ) ;
  }
  
  out = fcam_logprior_alpha(alpha, hyp_alpha1, hyp_alpha2) + Kplus * log(alpha) + lgamma(alpha) - 
            lgamma(N + alpha) + tmp ;
  return(out) ;
}

double fcam_MH_alpha(double current_par, double eps_alpha,
                     double hyp_alpha1, double hyp_alpha2,
                     int Kplus, int Kiter,
                     arma::vec cluster)
{
  double new_par ;
  double ratio ;
  new_par = exp(R::rnorm( log(current_par), eps_alpha )) ;
  
  ratio = exp( fcam_logpost_alpha(new_par, 
                                  hyp_alpha1, hyp_alpha2,
                                  cluster,
                                  Kplus, Kiter) - 
                fcam_logpost_alpha(current_par, 
                                   hyp_alpha1, hyp_alpha2,
                                   cluster,
                                   Kplus, Kiter) ) ;
  if(R::runif(0, 1) < ratio) { current_par = new_par ; }
  return(current_par) ;
}
