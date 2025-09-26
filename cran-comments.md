## SANple 0.2.0

In this release we

* Replaced `arma::is_finite` with `std::isfinite` to ensure compatibility with the latest `RcppArmadillo` release and CRAN policies.
* Implemented small code adjustments to address and resolve compiler warnings.
* Updated the MCMC functions to take a `burn in` value as input;
* Some algorithms ran `nrep-1` iterations. Updated to `nrep`;
* Improved efficiency of stick-breaking computation;
* Improved the initialization of the algorithms, streamlined some scripts;
* Changed `.cpp` `for loops` indexes from `int` to `unsigned int` when needed;
* Fixed a bug in the full conditional of the concentration parameter for the observational DP of CAM.


## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces no errors, warnings, or notes.

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
   
   New maintainer:
     Francesco Denti <francescodenti.personal@gmail.com>
   Old maintainer(s):
     Laura D'Angelo <laura.dangelo@live.com>

due to a change in the maintainer of the package.
      
- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... NOTE
Maintainer: 'Francesco Denti <francescodenti.personal@gmail.com>'

New maintainer:
  Francesco Denti <francescodenti.personal@gmail.com>
Old maintainer(s):
  Laura D'Angelo <laura.dangelo@live.com>
  
- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
