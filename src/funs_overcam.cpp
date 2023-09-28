#include "funs_overcam.h"

// log-prior of a vector of probabilities with Dirichlet prior distribution
double overcam_logprior_prob(arma::vec prob, double par)
{
  double out ;
  int K = prob.n_elem ;
  for(int i = 0; i < K; i++) { prob(i) = log(prob(i)) ; }
  out = lgamma( K * par ) - K * lgamma(par) + (par-1) * arma::accu(prob) ;
  return(out) ;
}

// log-prior on the hyperparameter of the Dirichlet distributions
// following Malsiner-Walli et al. (2016)
// alpha ~ Gamma(a, a*K) and beta ~ Gamma(b, b*L) 
double overcam_logprior_par_prob(double par, int dirichlet_dim, double hyperpar)
{
  double out = R::dgamma(par, hyperpar, 1/(hyperpar * dirichlet_dim), true) ;
  return(out) ;
}

// log-posterior distributional Dirichlet hyperparameter
double overcam_logpost_alpha(double alpha, arma::vec pi, double hyperpar)
{
  int K = pi.n_elem ;
  double out = overcam_logprior_prob(pi, alpha) + overcam_logprior_par_prob(alpha, K, hyperpar) ;
  return(out) ;
}

// log-posterior observational Dirichlet hyperparameter
double overcam_logpost_beta(double beta, arma::mat omega, double hyperpar)
{
  int L = omega.n_rows ;
  int K = omega.n_cols ;
  
  arma::vec tmp(K) ;
  for(int k = 0; k < K; k++) { tmp(k) = overcam_logprior_prob(omega.col(k), beta) ; }
  double out = arma::accu(tmp) + overcam_logprior_par_prob(beta, L, hyperpar) ;
  return(out) ;
}

// MH-step on distributional Dirichlet hyperparameter
double overcam_MH_alpha(double current_par, double eps, arma::vec pi, double hyperpar)
{
  double new_par ;
  double ratio ;
  new_par = R::rnorm( current_par, eps ) ;
  ratio = exp( overcam_logpost_alpha(new_par, pi, hyperpar) - overcam_logpost_alpha(current_par, pi, hyperpar) ) ;
  if(R::runif(0, 1) < ratio) current_par = new_par ;
  return(current_par) ;
}

// MH-step on observational Dirichlet hyperparameter
double overcam_MH_beta(double current_par, double eps, arma::mat omega, double hyperpar)
{
  double new_par ;
  double ratio ;
  new_par = R::rnorm( current_par, eps ) ;
  ratio = exp( overcam_logpost_beta(new_par, omega, hyperpar) - overcam_logpost_beta(current_par, omega, hyperpar) ) ;
  if(R::runif(0, 1) < ratio) current_par = new_par ;
  return(current_par) ;
}
