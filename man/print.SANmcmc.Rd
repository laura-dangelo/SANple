% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.SANmcmc}
\alias{print.SANmcmc}
\title{Print MCMC output}
\usage{
\method{print}{SANmcmc}(x, ...)
}
\arguments{
\item{x}{object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}},
\code{\link{sample_fSAN}}, or \code{\link{sample_CAM}}).}

\item{...}{ignored.}
}
\value{
The function prints a summary of the fitted model.
}
\description{
Print method for objects of class \code{SANmcmc}.
}
\seealso{
\code{\link{estimate_clusters}}, \code{\link{plot.SANmcmc}}
}
