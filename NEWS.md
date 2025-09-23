# SANple 0.2.0

* Major update of function names (C++ level) to uniform with accepted paper
* Removed fCAM model and its hybrid version (random concentration parameters)
* Replaced `arma::is_finite` with `std::isfinite` to ensure compatibility with the latest `RcppArmadillo` release and CRAN policies.

# SANple 0.1.1

* Updated the MCMC functions to take a `burn in` value as input;
* Some algorithms ran `nrep-1` iterations. Uniformed to `nrep`;
* Improved efficiency of stick-breaking computation;
* Improved the initialization of the algorithms, streamlined some scripts;
* Changed `.cpp` `for loops` indexes from `int` to `unsigned int` when needed;
* Fixed a bug in the full conditional of the concentration parameter for the observational DP of CAM.

# SANple 0.1.0

* Submitted to CRAN
* Released to the public the first official version
