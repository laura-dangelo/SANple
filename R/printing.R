
#' Print MCMC output
#' @description Print method for objects of class \code{SANmcmc}.
#'
#' @param x object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}},  
#' \code{\link{sample_fSAN}}, or \code{\link{sample_CAM}}).
#' @param ... ignored.
#'
#' @return The function prints a summary of the fitted model.
#' 
#' @seealso \code{\link{summary.SANmcmc}}, \code{\link{plot.SANmcmc}}
#' 
#' @export
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


#' Print cluster summary
#' @description Print the cluster-specific sample means and variances of the estimated observational and distributional partition.
#'
#' @param x object of class \code{SUMmcmc} (the result of a call to \code{\link{summary.SANmcmc}})
#' @param ... ignored.
#' 
#' @return The function prints a summary of the estimated clusters.
#' 
#' @export
print.summary_SANmcmc <- function(x, ...)
{
  
  n_dc <- length(unique(x$est_dc))
  
  cat("\n")
  cat("Summary of the estimated observational and distributional clusters \n\n")
  cat("----------------------------------\n")
  cat(paste("Estimated number of observational clusters:", length(unique(x$est_oc))),"\n")
  cat(paste("Estimated number of distributional clusters:", length(unique(x$est_dc))),"\n")
  
  cat("----------------------------------\n")
  
  for(j in 1:n_dc){
    cat(paste("\nDistributional cluster",j,"\n"))
    
    post_mean <- x$clus_means[,j]
    post_var <- x$clus_var[,j]
    
    occupied <- !is.na(x$clus_means[,j])
    Dsubj <- data.frame(cbind(post_mean, post_var))
    Dsubj <- Dsubj[occupied,]
    rownames(Dsubj) <- which(occupied)
    print(data.frame(Dsubj))
    
  }
}
