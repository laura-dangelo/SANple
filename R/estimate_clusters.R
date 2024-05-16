
#' Estimate observational and distributional clusters
#' @description Given the MCMC output, estimate the observational and distributional partitions using \code{\link[salso:salso]{salso::salso()}}.
#'
#' @param object object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, \code{\link{sample_fiSAN_sparsemix}}, 
#' \code{\link{sample_fSAN}}, \code{\link{sample_fSAN_sparsemix}}, or \code{\link{sample_CAM}}).
#' @param burnin the length of the burn-in to be discarded before estimating the clusters (default is 2/3 of the iterations).
#' @param ncores the number of CPU cores to use, i.e., the number of simultaneous runs at any given time. A value of zero indicates to use all cores on the system.
#'
#' @return Object of class \code{SANclusters}. The object contains:
#' @returns \code{est_oc} estimated partition at the observational level. It is an object of class \code{salso.estimate}.
#' @returns \code{est_dc} estimated partition at the distributional level. It is an object of class \code{salso.estimate}.
#' @returns \code{clus_means} cluster-specific sample means of the estimated partition.
#' @returns \code{clus_vars} cluster-specific sample variances of the estimated partition.
#' 
#' @seealso \code{\link[salso:salso]{salso::salso()}}, \code{\link{print.SANmcmc}}, \code{\link{plot.SANmcmc}}, \code{\link{print.SANclusters}}
#'
#' @examples 
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' out <- sample_fiSAN(nrep = 500, burn = 200,
#'                      y = y, group = g, 
#'                     nclus_start = 2,
#'                     maxK = 20, maxL = 20,
#'                     beta = 1)
#' estimate_clusters(out)
#'
#' @export
#' @importFrom salso salso
#' @useDynLib SANple
estimate_clusters <- function(object, burnin = 0, ncores = 0)
{
  
  if(burnin>0) { 
      OC <- object$sim$obs_cluster[-burnin,] 
      DC <- object$sim$distr_cluster[-burnin,]
  }else{
    OC <- object$sim$obs_cluster 
    DC <- object$sim$distr_cluster
  }
  
  estimated_oc <- suppressWarnings(salso::salso(OC, nCores = ncores)) 
  estimated_dc <- suppressWarnings(salso::salso(DC, nCores = ncores))
  
  n_oc <- length(unique(estimated_oc))
  n_dc <- length(unique(estimated_dc))
  
  means <- matrix(NA, n_oc, n_dc)
  vars <- matrix(0, n_oc, n_dc) # if singleton, variance set to zero
  
  dc_long <- estimated_dc[object$params$group]
  
  for(j in 1:n_dc) {
    
    suby <- object$params$y[dc_long == j]
    subcl <- estimated_oc[dc_long == j]
    
    for(k in unique(subcl)) {
      means[k,j] <- mean(suby[(subcl == k)]) 
      if(length(suby[subcl == k])>1) {
        vars[k,j] <- var(suby[subcl == k]) }
    }
  }
  

  
  out <- (list("est_oc" = estimated_oc,
              "est_dc" = estimated_dc,
              "clus_means" = means,
              "clus_vars" = vars))
  
  structure(out, class = c("SANclusters", class(out)))
}
