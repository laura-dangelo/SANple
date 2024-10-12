#' #' Plotting MCMC output
#' #' @description Plot method for objects of class \code{SANmcmc}. 
#' #' The function displays two graphs, meant to provide a quick assessment of the estimated distributional clusters.
#' #' 
#' #' @param x object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, 
#' #' \code{\link{sample_fSAN}}, or \code{\link{sample_CAM}}).
#' #' @param estimated_clusters the output of a call to \code{\link{estimate_clustering_mcmc}} (optional). It can be used to speed up the function if the partition has already been computed. 
#' #' If \code{estimated_clusters = NULL}, the displayed partition is computed using \code{\link{estimate_clustering_mcmc}}.
#' #' @param burnin the length of the burn-in to be discarded (default is 2/3 of the iterations).
#' #' @param palette_brewed (logical) the color palette to be used. Default is \code{R} base colors (\code{palette_brewed = FALSE}).
#' #' @param ncores if the partition is computed, the number of CPU cores to use to estimate the clusters, i.e., the number of simultaneous runs at any given time. A value of zero indicates to use all cores on the system.
#' #' @param ... additional graphical parameters to be passed when \code{type = "scatter"} is used.
#' #'
#' #' @return The function plots a summary of the fitted model.
#' #'
#' #' @seealso \code{\link{print.SANmcmc}}, \code{\link{estimate_clustering_mcmc}}
#' #'
#' #' @examples 
#' #' set.seed(123)
#' #' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' #' g <- c(rep(1,30), rep(2, 30))
#' #' out <- SANple::sample_fiSAN(y = y, group = g)
#' #' plot(out, type = "ecdf", palette_brewed = TRUE)
#' #'
#' #' @importFrom graphics abline lines points boxplot par
#' #' @importFrom grDevices colorRampPalette
#' #' @importFrom scales alpha
#' #' @importFrom stats ecdf
#' #' @export
#' #' 
#' #' Voglio mettere solo DC in plot di SANmcmc, come per SANvi
#' plot.SANmcmc <- function(x,
#'                          params, 
#'                          show_density = TRUE, 
#'                          show_burnin = TRUE, 
#'                          length_burnin = NULL, 
#'                          show_convergence = TRUE, 
#'                          trunc_plot = 10, ...) 
#' {
#'   oldpar <- par(no.readonly = TRUE) 
#'   on.exit(par(oldpar)) 
#'   
#'   if(show_density){
#'     .traces_and_density(object$sim, params, 
#'                         show_burnin,
#'                         length_burnin,
#'                         show_convergence,
#'                         trunc_plot)
#'   } else {
#'     .traces(object$sim, params, 
#'             show_burnin,
#'             length_burnin,
#'             show_convergence,
#'             trunc_plot)
#'   }
#'   devAskNewPage(ask = F)
#'   
#' }
#' 
