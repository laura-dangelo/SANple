
#' Print MCMC output
#' @description Print method for objects of class \code{SANmcmc}.
#'
#' @param x object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, \code{\link{sample_fiSAN_sparsemix}}, 
#' \code{\link{sample_fSAN}}, \code{\link{sample_fSAN_sparsemix}}, or \code{\link{sample_CAM}}).
#' @param ... ignored.
#'
#' @return The function prints a summary of the fitted model.
#' 
#' @seealso \code{\link{estimate_clusters}}, \code{\link{plot.SANmcmc}}
#' 
#' @export
#' @useDynLib SANple
print.SANmcmc <- function(x, ...)
{
  cat("\n")
  cat(paste("MCMC result of", x$model, "model \n"))
  cat("-----------------------------------------------\n")
  cat(paste("Model estimated on", length(x$params$y), "total observations and",  length(unique(x$params$group)), "groups \n"))
  cat(paste("Total MCMC iterations:", x$params$nrep, "\n"))
  cat(paste("maxL:",x$params$maxL,"- maxK:",x$params$maxK,"\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time[[1]]),3), attr(x$time, "units"),"\n\n"))
}