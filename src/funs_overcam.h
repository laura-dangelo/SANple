#ifndef OVERCAMFUNS_H
#define OVERCAMFUNS_H

#include "funs_san.h"

double overcam_logprior_prob(arma::vec prob, double par) ;

double overcam_logprior_par_prob(double par, int dirichlet_dim, double hyperpar) ;

double overcam_logpost_alpha(double alpha, arma::vec pi, double hyperpar) ;

double overcam_logpost_beta(double beta, arma::mat omega, double hyperpar) ;

double overcam_MH_alpha(double current_par, double eps, arma::vec pi, double hyperpar) ;

double overcam_MH_beta(double current_par, double eps, arma::mat omega, double hyperpar) ;

#endif

