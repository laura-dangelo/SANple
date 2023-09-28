#ifndef FICAMFUN_H
#define FICAMFUN_H

#include "funs_cam.h"
#include "funs_fcam.h"
#include "funs_san.h"

Rcpp::List ficam_weights_update_slice_sampler(const arma::vec& group, 
                                              arma::vec S_iter, 
                                              double alpha,
                                              arma::vec xi,
                                              int maxK) ;
#endif

