#' Summary of model: estimate observational and distributional clusters
#' 
#' @description Given the MCMC output, estimate the observational and distributional partitions using \code{\link[salso:salso]{salso::salso()}}.
#'
#' @param object object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, 
#' \code{\link{sample_fSAN}}, or \code{\link{sample_CAM}}).
#' @param burnin the length of the burn-in to be discarded before estimating the clusters (default is 2/3 of the iterations).
#' @param ncores the number of CPU cores to use, i.e., the number of simultaneous runs at any given time. A value of zero indicates to use all cores on the system.
#'
#' @return a list of class \code{clustering} containing 
#' \itemize{
#'   \item \code{obs_level}: a data frame containing the data values, their group indexes, the observational and distributional clustering assignments for each observation.
#'   \item \code{dis_level}: a vector with the distributional clustering assignment for each unit.
#' }
#'  
#' @seealso \code{\link[salso:salso]{salso::salso()}}, \code{\link{print.SANmcmc}}, \code{\link{plot.SANmcmc}}
#'
#' @examples 
#' set.seed(1234)
#' y <- c(rnorm(400,0,0.3), rnorm(200,5,0.3))
#' g <- c(rep(1,300), rep(2, 300))
#' object <- sample_fiSAN(y = y, group = g)
#' summary(object)
#'
#' @export
#' @importFrom salso salso
summary.SANmcmc <- function(object, burnin = 0, ncores = 1, ...)
{
  
  if(burnin>0) { 
      OC <- object$sim$obs_cluster[-burnin,] 
      DC <- object$sim$distr_cluster[-burnin,]
  }else{
    OC <- object$sim$obs_cluster 
    DC <- object$sim$distr_cluster
  }
  estimated_oc <- suppressWarnings(salso::salso(OC, nCores = ncores,maxNClusters = 
                                                  max(object$params$group))) 
  estimated_dc <- suppressWarnings(salso::salso(DC, nCores = ncores,maxNClusters = 
                                                  max(object$params$group)))

  D <- data.frame(Y  = object$params$y,
                  G  = object$params$group,
                  OC = estimated_oc,
                  DC = rep(estimated_dc, 
                           object$params$Nj))
  
  L <- (list(obs_level = D, 
             dis_level = estimated_dc))
  
  structure(L,class = c("SUMmcmc",class(L)))
  
}


#' @name summary.SANmcmc
#' 
#' @param x an object of class \code{mcmc_clustering}, which can be obtained from the
#' function \code{\link{summary.SANmcmc}}.
#' @param DC_num an integer or a vector of integers indicating which distributional clusters to plot.
#' @param type what type of plot should be drawn (only for the left-side plot). Possible types are "boxplot", "ecdf", and "scatter". 
#' @param palette_brewed (logical) the color palette to be used. Default is \code{R} base colors (\code{palette_brewed = FALSE}).
#' @param ... ignored.
#'
#' @importFrom graphics abline lines points boxplot
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
plot.SUMmcmc <- function(x,
                         DC_num = NULL,
                         type = c("ecdf", "boxplot", "scatter"),
                         palette_brewed = FALSE,
                         ...) {
  
  type   <- match.arg(type)
  if (is.null(DC_num)) {
    DC_num <- 1:max(x$dis_level)
  }
  cat(DC_num)
  if(max(DC_num) > max(x$dis_level)){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", max(x$dis_level) ))
  }
  
  max_CD <- max(x$dis_level)
  if(palette_brewed){
    colpal <-
      rev(
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(max_CD))
  }else{
    colpal <- 1:max_CD
  }
  
  dix <- rank(tapply(x$obs_level$Y, x$obs_level$DC, stats::median))
  
  ind_ord_dis <- sapply(x$dis_level, function(g) dix[g])
  ind_ord_obs <- dix[x$obs_level$DC]
  
  
  inds_row <- which(ind_ord_obs %in% DC_num)
  inds_col <- which(ind_ord_dis %in% DC_num)
  
  suby <- (x$obs_level$Y[inds_row])
  subg <- (x$obs_level$G[inds_row])
  subDC <- x$obs_level$DC[inds_row]
  subOC <- x$obs_level$OC[inds_row]
  
  if(length(DC_num) == max(x$dis_level)){
    main_title <- paste0("All distributional clusters")
  }else if( length(DC_num) > 8 ){
    main_title <- NULL
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }
  
  
  
  if (type == "ecdf") {
    X <- unique(x$obs_level$G[inds_row])
    Nj <- table(x$obs_level$G)
    ysteps <- (0:(Nj[inds_col][1] - 1)) / 
      (Nj[inds_col][1] - 1)
    xsteps <- sort(suby[subg == X[1]])
    
    
    plot(
      ysteps ~ xsteps,
      type = "l",
      xlim = c(min(suby), max(suby)),
      ylim = c(0, 1),
      col = scales::alpha(colpal[ind_ord_dis[inds_col][1]], .5),
      xlab = "y",
      ylab = "eCDF",
      main = paste0("eCDFs colored by DC\n",main_title)
    )
    graphics::points(
      (ysteps ~ xsteps),
      cex = .1,
      col = scales::alpha(colpal[ind_ord_dis[inds_col][1]], .5)
    )
    graphics::abline(h = c(0, 1),
                     col = "gray",
                     lty = 3)
    
    
    for (j in 2:length(X)) {
      ysteps = (0:(Nj[X[j]] - 1)) / (Nj[X[j]] - 1)
      xsteps = sort(suby[subg == X[j]])
      
      graphics::lines((ysteps ~ xsteps),
                      col = scales::alpha(colpal[ind_ord_dis[inds_col][j]], .5))
      graphics::points(
        (ysteps ~ xsteps),
        cex = .1,
        col = scales::alpha(colpal[ind_ord_dis[inds_col][j]], .5)
      )
      
    }
    
  } else if (type=="boxplot"){
    
    graphics::boxplot(
      suby ~ subg,
      col = scales::alpha(colpal[ind_ord_dis[inds_col]], .7),
      pch = ".",
      main = paste0("Boxplots colored by DC\n",main_title),
      ylab = "y",
      xlab = "group"
    )
  }else{
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par))    
    
    graphics::par(mfrow=c(1,2))
    plot(suby ~ jitter(subg),
         pch=".",
         col=subDC,
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by DC\n",main_title)
    )
    plot(suby ~ jitter(subg),
         pch=".",
         col=subOC,
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by OC\n",main_title)
    )
  }
}

#' @name summary.SANmcmc
#'
#' @export
#'
print.SUMmcmc <- function(x, ...){
  cat(paste("Number of estimated OCs:", length(unique(x$obs_level$OC))),"\n")
  cat(paste("Number of estimated DCs:",length(unique(x$dis_level))))
  invisible(x)
}




