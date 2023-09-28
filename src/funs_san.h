#ifndef SANFUNCTIONS_H
#define SANFUNCTIONS_H

#include "common_functions.h"

arma::vec dirichlet_sample_distr_weights(arma::vec S_iter, 
                                         double alpha, 
                                         int K_iter, 
                                         int maxK) ;

arma::mat dirichlet_sample_obs_weights(arma::vec M_iter, 
                                       arma::vec clusterD_long,
                                       double beta,
                                       int K_iter, int L_iter,
                                       int maxK, int maxL) ;

arma::vec sample_distr_cluster(const arma::vec& y, const arma::vec& group,
                               arma::vec pi, arma::mat omega,
                               int K_iter, int L_iter,
                               arma::vec mu, arma::vec sigma2) ;

arma::vec sample_obs_cluster(const arma::vec& y, 
                             arma::vec clusterD_long,
                             arma::mat omega,
                             int K_iter, int L_iter,
                             arma::vec mu, arma::vec sigma2); 

#endif

