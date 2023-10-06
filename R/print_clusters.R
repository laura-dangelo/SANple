#' Print cluster summary
#' @description Print the cluster-specific sample means and variances of the estimated observational and distributional partition.
#'
#' @param x object of class \code{SANclusters} (the result of a call to \code{\link{estimate_clusters}})
#' @param ... ignored.
#' 
#' @return The function prints a summary of the estimated clusters.
#' 
#' @export
print.SANclusters <- function(x, ...)
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