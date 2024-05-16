## SANple 0.1.1

In this release we

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
  Maintainer: ‘Laura D'Angelo <laura.dangelo@live.com>’
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1080/01621459.2021.1933499
      From: DESCRIPTION
      Status: Forbidden
      Message: 403
      
However, the links work when building the documentation.

- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 0 note

- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
