
#' Traceplot: plot MCMC chains
#' @description Check the convergence of the MCMC through visual inspection of the chains.
#' 
#' @usage 
#' traceplot(object, params, 
#'           show_density = TRUE, 
#'           show_burnin = TRUE, 
#'           length_burnin = NULL, 
#'           show_convergence = TRUE, 
#'           trunc_plot = 10)
#'
#' @param object object of class \code{SANmcmc} (the result of a call to \code{\link{sample_fiSAN}}, 
#' \code{\link{sample_fiSAN_sparsemix}}, \code{\link{sample_fSAN}}, \code{\link{sample_fSAN_sparsemix}}, or \code{\link{sample_CAM}}).
#' @param params vector of strings with the names of the parameters to check.
#' @param show_burnin logical (default \code{TRUE}). Whether the first part of the chains should be plotted in the traceplots.
#' @param show_density  logical (default \code{TRUE}). Whether a kernel estimate of the density should be plotted. The burn-in is always discarded.
#' @param length_burnin if \code{show_burnin = FALSE}, the length of the burn-in to be discarded.
#' @param show_convergence logical (default \code{TRUE}). Whether a superimposed red line of the cumulative mean should be plotted.
#' @param trunc_plot integer (default = 10). For multidimensional parameters, the maximum number of components to be plotted.
#' 
#' @note The function is not available for the observational weights \eqn{\omega}.
#' 
#' @return The function displays the traceplots of the MCMC algorithm.
#' 
#' @examples 
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' out <- sample_fiSAN(nrep = 500, y = y, group = g, 
#'                     nclus_start = 2,
#'                     maxK = 20, maxL = 20,
#'                     beta = 1)
#' traceplot(out, params = c("mu", "sigma2"), trunc_plot = 2)
#' 
#' 
#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @export
#' @useDynLib SANple
traceplot <- function(object, params, 
                     show_density = TRUE,
                     show_burnin = TRUE,
                     length_burnin = NULL,
                     show_convergence = TRUE,
                     trunc_plot = 10)
{
  
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar)) 

  if(show_density){
    .traces_and_density(object$sim, params, 
            show_burnin,
            length_burnin,
            show_convergence,
            trunc_plot)
  } else {
    .traces(object$sim, params, 
            show_burnin,
            length_burnin,
            show_convergence,
            trunc_plot)
  }
  devAskNewPage(ask = F)
}



#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @keywords internal
.traces <- function(sim, params,
                    show_burnin,
                    length_burnin,
                    show_convergence,
                    trunc_plot)
{
  burnin <- 1
  if((!show_burnin) & (!is.null(length_burnin))) {burnin <- 1:length_burnin}
  if((!show_burnin) & (is.null(length_burnin))) {burnin <- 1:floor(nrow(sim$mu)/3)}
  if("omega" %in% params) {
    warning("Traceplot for omega is not available")
    params <- params[params!="omega"]
  }

  stringg <- paste0("sim$", params)
  count <- 0
  tmp_names <- c()
  
  for(i in 1:length(params))
  {
    tmp <- eval(parse(text=stringg[i]))
    if(dim(tmp)[2] > trunc_plot) { tmp_names <- c(tmp_names, params[i]) }
    for(j in 1:min(dim(tmp)[2],  trunc_plot)){
      count <- count+1
    }
  }
  par(mfrow = c(min(3,count),1))

  devAskNewPage(ask <- F)
  for(i in 1:length(params))
  {
    tmp <- eval(parse(text=stringg[i]))
    for(j in 1:min(dim(tmp)[2],  trunc_plot)){
      
      if(dim(tmp)[2] == 1) { namem <- paste0(params[i]) } else { namem <- paste0(params[i],"_",j) }
      
      plot(tmp[-burnin,j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
      if(show_convergence) {
        lines(1:length(tmp[-burnin,j]), cumsum((tmp[-burnin,j]))/(1:length(tmp[-burnin,j])), col=2)
      }
      devAskNewPage(ask = T)
    }
  }
  print(paste0("Output truncated at ", trunc_plot, " for ", paste0(tmp_names, collapse = ", "), "."))
}




#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @importFrom stats density
#' @keywords internal
.traces_and_density <- function(sim, params,
                   show_burnin,
                   length_burnin,
                   show_convergence,
                   trunc_plot)
{
  burnin <- 1
  if((!show_burnin) & (!is.null(length_burnin))) {burnin <- 1:length_burnin}
  if((!show_burnin) & (is.null(length_burnin))) {burnin <- 1:floor(nrow(sim$mu)/3)}
  if("omega" %in% params) {
    warning("Traceplot for omega is not available")
    params <- params[params!="omega"]
  }
  
  stringg <- paste0("sim$", params)
  count <- 0
  tmp_names <- c()
  
  for(i in 1:length(params))
  {
    tmp <- eval(parse(text=stringg[i]))
    if(dim(tmp)[2] > trunc_plot) { tmp_names <- c(tmp_names, params[i]) }
    for(j in 1:min(dim(tmp)[2],  trunc_plot)){
      count <- count+1
    }
  }
  par(mfrow = c(min(3,count),2))
  
  devAskNewPage(ask = F)
  for(i in 1:length(params))
  {
    tmp <- eval(parse(text=stringg[i]))
    for(j in 1:min(dim(tmp)[2],  trunc_plot)){
      
      if(dim(tmp)[2] == 1) { namem <- paste0(params[i]) } else { namem <- paste0(params[i],"_",j) }
      
      plot(tmp[-burnin,j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
      if(show_convergence) {
        lines(1:length(tmp[-burnin,j]), cumsum((tmp[-burnin,j]))/(1:length(tmp[-burnin,j])), col=2)
      }
      
      plot(density(tmp[-c(1:max(length_burnin, floor(nrow(sim$mu)/3) )),j]), main = namem)
      
      devAskNewPage(ask = T)
    }
  }
  print(paste0("Output truncated at ", trunc_plot, " for ", paste0(tmp_names, collapse = ", "), "."))
}

