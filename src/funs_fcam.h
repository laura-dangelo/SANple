#ifndef FCAMFUN_H
#define FCAMFUN_H

#include "funs_san.h"

arma::vec fcam_sample_distr_cluster(const arma::vec& y, const arma::vec& group,
                                    arma::vec pi, arma::mat omega,
                                    arma::vec M_iter,
                                    int K_iter) ;

arma::vec fcam_sample_obs_cluster(const arma::vec& y, 
                                   arma::vec clusterD_long,
                                   arma::mat omega,
                                   int K_iter, int L_iter,
                                   arma::vec mu, arma::vec sigma2) ;

double lbeta(double a, double b) ;

double logprior_maxx(int maxx, double hyp_max1,
                     double hyp_max2, double hyp_max3) ;

double logpost_maxx(int Kiter, 
                    double hyp_max1, 
                    double hyp_max2, 
                    double hyp_max3, 
                    double alpha,
                    arma::vec& cluster,
                    int Kplus) ;

int fcam_sample_K(int maxK,
                  int Kplus,
                  double hyp_max1, 
                  double hyp_max2, 
                  double hyp_max3,
                  double d_par,
                  arma::vec cluster) ;

double fcam_logprior_alpha(double alpha, 
                           double hyp_alpha1, double hyp_alpha2) ;

double fcam_logpost_alpha(double alpha, 
                          double hyp_alpha1, double hyp_alpha2,
                          arma::vec cluster,
                          int Kplus,
                          int Kiter) ;

double fcam_MH_alpha(double current_par, double eps_alpha,
                     double hyp_alpha1, double hyp_alpha2,
                     int Kplus, int Kiter,
                     arma::vec cluster) ;
  
#endif

