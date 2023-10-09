Please note this is a resubmission of this package to CRAN.

Thank you for checking our package. After the inspection, we were notified of the following:

* Please do not start the description with "This package", package name, title or similar.

The new incipit in the description is: "Estimate Bayesian nested mixture models via Markov Chain Monte Carlo methods.[...]".

* Please add some form of linking to your references.

We added the links for the papers we mentioned.

* Please always explain all acronyms in the description text. -> MCMC

Done.

* Please call the on.exit() functions right after defining 'oldpar' and change the option only afterwards.

We corrected the code as requested.

* Please unwrap the examples if that is feasible and if they can be 
executed in < 5 sec for each Rd file or create additionally small toy 
examples to allow automatic testing.

All the old examples have been reduced and are now testable. On our machines, the tests took less than a second to run.

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces no errors, warnings, or notes.

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Laura D'Angelo <laura.dangelo@live.com>’
  
  New submission

and

Found the following (possibly) invalid DOIs:
    DOI: 10.1080/01621459.2021.1933499
      From: DESCRIPTION
      Status: Forbidden
      Message: 403
    DOI: 10.1111/biom.13626
      From: DESCRIPTION
      Status: Forbidden
      Message: 403

However, the links work when building the documentation.

- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Laura D'Angelo <laura.dangelo@live.com>'

  New submission

- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
