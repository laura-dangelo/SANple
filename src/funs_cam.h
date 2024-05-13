#ifndef CAMFUN_H
#define CAMFUN_H

#include "common_functions.h"

double fun_xi(double kappa, int i) ;

int compute_trunc(double u_min, double kappa) ;
  
arma::vec stick_breaking(arma::vec beta_var) ;

Rcpp::List weights_update_slice_sampler(const arma::vec& y, const arma::vec& group, 
                                        arma::vec S_iter, arma::vec M_iter, 
                                        arma::vec clusterD_long,
                                        arma::vec xi,
                                        double & alpha, double & beta, 
                                        int & maxK, int & maxL) ;

arma::vec slicedDP_sample_distr_cluster(const arma::vec& group,
                                        arma::vec M_iter,
                                        arma::vec pi, arma::mat omega,
                                        arma::vec u_D,
                                        arma::vec xi,
                                        int maxK_iter) ;
// This function performs step 5 of Algorithm 1 of Denti et al. (2021):
// sample distributional cluster allocation S_j
arma::vec slicedDP_sample_distr_cluster2(const arma::vec& group,
                                         arma::vec M_iter,
                                         arma::vec pi, arma::mat omega,
                                         arma::vec u_D,
                                         arma::vec xi,
                                         int maxK_iter);

arma::vec slicedDP_sample_obs_cluster(const arma::vec& y, 
                                      arma::vec clusterD_long,
                                      arma::mat omega,
                                      arma::vec u_O,
                                      arma::vec xi,
                                      int maxL_iter,
                                      arma::vec mu, arma::vec sigma2) ;

double sample_alpha(double old_alpha,
                    double hyp_alpha1, double hyp_alpha2,
                    int G , arma::vec S_iter) ;
  
double sample_beta(double hyp_beta1, double hyp_beta2,
                   arma::vec S_iter, arma::vec M_iter, 
                   arma::vec clusterD_long,
                   arma::mat obs_beta_rv) ;
#endif

