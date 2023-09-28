#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <RcppArmadillo.h>

arma::vec rdirichlet(arma::vec par) ;

int sample_i(arma::vec ids, arma::vec prob) ;

arma::vec relabel_arma(arma::vec cluster) ;

Rcpp::List sample_model_parameters(const arma::vec& y, 
                                   arma::vec M_iter,
                                   int maxL,
                                   double m0, double tau0,
                                   double lambda0, double gamma0) ;

#endif

