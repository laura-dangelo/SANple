
#' Plotting MCMC output
#' @description Plot method for objects of class \code{SANmcmc}. 
#' The function displays two graphs, meant to analyze the estimated distributional and observational clusters.
#' 
#' @param x object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, \code{\link{sample_fiSAN_sparsemix}}, 
#' \code{\link{sample_fSAN}}, \code{\link{sample_fSAN_sparsemix}}, or \code{\link{sample_CAM}}).
#' @param type what type of plot should be drawn (only for the left-side plot). Possible types are "boxplot", "ecdf", and "scatter". 
#' @param estimated_clusters the output of a call to \code{\link{estimate_clusters}} (optional). It can be used to speed up the function if the partition has already been computed. 
#' If \code{estimated_clusters = NULL}, the displayed partition is computed using \code{\link{estimate_clusters}}.
#' @param burnin the length of the burn-in to be discarded (default is 2/3 of the iterations).
#' @param palette_brewed (logical) the color palette to be used. Default is \code{R} base colors (\code{palette_brewed = FALSE}).
#' @param ncores if the partition is computed, the number of CPU cores to use to estimate the clusters, i.e., the number of simultaneous runs at any given time. A value of zero indicates to use all cores on the system.
#' @param ... additional graphical parameters to be passed when \code{type = "scatter"} is used.
#'
#' @return The function plots a summary of the fitted model.
#'
#' @seealso \code{\link{print.SANmcmc}}, \code{\link{estimate_clusters}}
#'
#' @examples 
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' out <- sample_fiSAN(nrep = 500, y = y, group = g, 
#'                     nclus_start = 2,
#'                     maxK = 20, maxL = 20,
#'                     beta = 1)
#' plot(out, type = "ecdf", palette_brewed = TRUE)
#'
#' @importFrom graphics abline lines points boxplot par
#' @importFrom grDevices colorRampPalette
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats ecdf
#' @export
plot.SANmcmc <- function(x,
                         type = c("boxplot", "ecdf", "scatter"),
                         estimated_clusters = NULL,
                         burnin = NULL,
                         palette_brewed = FALSE, ncores = 0, ...) 
{
  type <- match.arg(type)
  
  if(!is.null(estimated_clusters)) { 
    burnin <- 1:(x$params$nrep - nrow(attr(estimated_clusters$est_oc, "draws")))
    estimated_oc <- estimated_clusters$est_oc
    estimated_dc <- estimated_clusters$est_dc
  } else {
    if(is.null(burnin)) { burnin <- 1:round(x$params$nrep/3*2) 
    } else { burnin <- 1:burnin }
  }
  
  if(is.null(estimated_clusters)) { 
    estimated_oc <- suppressWarnings(salso::salso(x$sim$obs_cluster[-burnin,], nCores = ncores)) 
    estimated_dc <- suppressWarnings(salso::salso(x$sim$distr_cluster[-burnin,], nCores = ncores)) 
  }
  
  posterior_means <- tapply(x$params$y, estimated_oc, mean)
  
  
  max_CD <- max(estimated_dc)
  max_OC <- max(estimated_oc)
  if(palette_brewed){
    colpal <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(max_CD))
    colpal2 <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(max_OC))
  }else{
    colpal <- 1:max_CD
    colpal2 <- 1:max_OC
  }
  
  
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar)) 
  
  par(mfrow=c(1,2))
  
  if(type == "ecdf") {  
    
    ecdfs <- list()
    DCs <- c()
    for(j in 1:length(unique(x$params$group))) {
      idj <- unique(x$params$group)[j]
      ecdfs[[j]] <- ecdf(x$params$y[x$params$group == idj])
      DCs[j] <- estimated_dc[idj]
    }
    plot(ecdfs[[1]], verticals=TRUE, do.points=FALSE, col = colpal[DCs[1]], 
         xlim = range(x$params$y), main = "eCDFs colored by DC",
         xlab = "y", ylab = "eCDF",)
    for(j in 2:length(ecdfs)) { plot(ecdfs[[j]], verticals=TRUE, do.points=FALSE, add=TRUE, col = colpal[DCs[j]]) }
    
  } else if(type=="boxplot") {
    
    graphics::boxplot(
      x$params$y ~ x$params$group,
      col = scales::alpha(colpal[estimated_dc], .7),
      main = paste0("Boxplots colored by DC"),
      ylab = "y",
      xlab = "Group"
    )
    
  } else if(type=="scatter") {
    
    graphics::par(mfrow=c(1,2))
    plot(x$params$y ~ jitter(x$params$group),
         col=colpal[estimated_dc[x$params$group]],
         xlab = "Group index",
         ylab = "y",
         main = "Observations colored by DC", ...)
  }
  
  plot(x$params$y ~ jitter(x$params$group),
       # pch=".",
       col=colpal2[estimated_oc],
       xlab = "Group index",
       ylab = "y",
       main = "Observations colored by OC", ...)
  abline(h = posterior_means, col=4, lty=2)
  
}

