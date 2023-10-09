#' Sample fiSAN
#' 
#' @description \code{sample_fiSAN} is used to perform posterior inference under the finite-infinite shared atoms nested (fiSAN) model with Gaussian likelihood. 
#' The model uses a Dirichlet process mixture prior at the distributional level, 
#' and Dirichlet mixture with an unknown number of components (following the specification of Frühwirth-Schnatter et al., 2021) at the observational level.
#' The algorithm for the nonparametric component is based on the slice sampler for DPM of Kalli, Griffin and Walker (2011).
#' 
#' 
#' @usage 
#' sample_fiSAN(nrep, y, group, 
#'              maxK = 50, maxL = 50, 
#'              m0 = 0, tau0 = 0.1, lambda0 = 3, gamma0 = 2,
#'              hyp_alpha1 = 1, hyp_alpha2 = 1,
#'              hyp_beta1 = 6, hyp_beta2 = 3, 
#'              eps_beta = NULL,
#'              alpha = NULL, beta = NULL,
#'              warmstart = TRUE, nclus_start = NULL,
#'              mu_start = NULL, sigma2_start = NULL,
#'              M_start = NULL, S_start = NULL,
#'              alpha_start = NULL, beta_start = NULL,
#'              progress = TRUE, seed = NULL)
#' 
#' @param nrep Number of MCMC iterations.
#' @param y Vector of observations.
#' @param group Vector of the same length of y indicating the group membership (numeric).
#' @param maxK Maximum number of distributional clusters \eqn{K} (default = 50).
#' @param maxL Maximum number of observational clusters \eqn{L} (default = 50).
#' @param m0,tau0,lambda0,gamma0 Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}. Default is (0, 0.1, 3, 2).
#' @param hyp_alpha1,hyp_alpha2 If a random \eqn{\alpha} is used, (\code{hyp_alpha1},\code{hyp_alpha2}) specify the hyperparameters (default = (1,1)).
#' The prior is \eqn{\alpha} ~ Gamma(\code{hyp_alpha1}, \code{hyp_alpha2}).
#' @param hyp_beta1,hyp_beta2,eps_beta If a random \eqn{\beta} is used, (\code{hyp_beta1},\code{hyp_beta2}) specifies the hyperparameter (default = (6,3)). 
#' The prior is \eqn{\beta\sim Gamma}(\code{hyp_beta1}, \code{hyp_beta2}). In this case, \code{eps_beta} is the tuning parameter of the MH step.
#' @param alpha Distributional DP parameter if fixed (optional). The distribution is \eqn{\pi\sim GEM (\alpha)}.
#' @param beta Observational Dirichlet parameter if fixed (optional). The distribution is Dirichlet( \code{rep(beta/maxL, maxL)} ).
#' @param warmstart,nclus_start Initialization of the observational clustering. 
#' \code{warmstart} is logical parameter (default = \code{TRUE}) of whether a kmeans clustering should be used to initialize the chains.
#' An initial guess of the number of observational clusters can be passed via the \code{nclus_start} parameter (optional)
#' @param mu_start,sigma2_start,M_start,S_start,alpha_start,beta_start Starting points of the MCMC chains (optional). Default is \code{nclus_start = min(c(maxL, 30))}.
#' \code{mu_start, sigma2_start} are vectors of length \code{maxL}. 
#' \code{M_start} is a vector of observational cluster allocation of length N.
#' \code{S_start} is a vector of observational cluster allocation of length J.
#' \code{alpha_start, alpha_start} are numeric. 
#' @param progress show a progress bar? (logical, default TRUE).
#' @param seed set a fixed seed.
#'
#'
#' @details 
#' \strong{Data structure}
#' 
#' The finite-infinite common atoms mixture model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups. 
#' The data should be continuous observations \eqn{(Y_1,\dots,Y_J)}, where each \eqn{Y_j = (y_{1,j},\dots,y_{n_j,j})} 
#' contains the \eqn{n_j} observations from group \eqn{j}, for \eqn{j=1,\dots,J}. 
#' The function takes as input the data as a numeric vector \code{y} in this concatenated form. Hence \code{y} should be a vector of length
#' \eqn{n_1+\dots+n_J}. The \code{group} parameter is a numeric vector of the same size as \code{y} indicating the group membership for each
#' individual observation. 
#' Notice that with this specification the observations in the same group need not be contiguous as long as the correspondence between the variables
#' \code{y} and \code{group} is maintained.
#'
#' \strong{Model}
#' 
#' The data are modeled using a univariate Gaussian likelihood, where both the mean and the variance are observational-cluster-specific, i.e., 
#' \deqn{y_{i,j}\mid M_{i,j} = l \sim N(\mu_l,\sigma^2_l)}
#' where \eqn{M_{i,j} \in \{1,\dots,L \}} is the observational cluster indicator of observation \eqn{i} in group \eqn{j}.
#' The prior on the model parameters is a Normal-Inverse-Gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)}, 
#' i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim Gamma(\lambda_0, \gamma_0)} (shape, rate).
#'
#' \strong{Clustering}
#' 
#' The model performs a clustering of both observations and groups. 
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,2,\dots\}}, with 
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,2,\dots}
#' The distribution of the probabilities is \eqn{ \{\pi_k\}_{k=1}^{\infty} \sim GEM(\alpha)}, 
#' where GEM is the Griffiths-Engen-McCloskey distribution of parameter \eqn{\alpha}, 
#' which characterizes the stick-breaking construction of the DP (Sethuraman, 1994). 
#' 
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,2,\dots \, ; \: l = 1,\dots,L. }
#' The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim Dirichlet_L(\beta/L,\dots,\beta/L)} for all \eqn{k = 1,2,\dots}. 
#' Moreover, the dimension \eqn{L} is random (see Frühwirth-Schnatter et al., 2021).
#'
#'
#' @return \code{sample_fiSAN} returns four objects:
#' \itemize{
#'   \item \code{model}: name of the fitted model.
#'   \item \code{params}: list containing the data and the parameters used in the simulation. Details below.
#'   \item \code{sim}: list containing the simulated values (MCMC chains). Details below.
#'   \item \code{time}: total computation time.
#' }
#' 
#'  
#' \strong{Data and parameters}:
#' \code{params} is a list with the following components:
#' \describe{
#' \item{\code{nrep}}{Number of MCMC iterations.}
#' \item{\code{y, group}}{Data and group vectors.}
#' \item{\code{maxK, maxL}}{Maximum number of distributional and observational clusters.}
#' \item{\code{m0, tau0, lambda0, gamma0}}{Model hyperparameters.}
#' \item{(\code{hyp_alpha1,hyp_alpha2}) or \code{alpha}}{Either the hyperparameters on \eqn{\alpha} (if \eqn{\alpha} random), or the value for \eqn{\alpha} (if fixed).}
#' \item{(\code{hyp_beta1, hyp_beta2, eps_beta}) or \code{beta}}{Either the hyperparameters on \eqn{\beta} and MH step size (if \eqn{\beta} random), or the value for \eqn{\beta} (if fixed).}
#' }
#' 
#' 
#' \strong{Simulated values}:
#' \code{sim} is a list with the following components:
#' \describe{
#' \item{\code{mu}}{Matrix of size (\code{nrep}, \code{maxL}).
#'    Each row is a posterior sample of the mean parameter for each observational cluster \eqn{(\mu_1,\dots\mu_L)}.}
#' \item{\code{sigma2}}{Matrix of size (\code{nrep}, \code{maxL}). 
#'     Each row is a posterior sample of the variance parameter for each observational cluster \eqn{(\sigma^2_1,\dots\sigma^2_L)}.}
#' \item{\code{obs_cluster}}{Matrix of size (\code{nrep}, n), with n = \code{length(y)}. 
#'    Each row is a posterior sample of the observational cluster allocation variables \eqn{(M_{1,1},\dots,M_{n_J,J})}. }
#' \item{\code{distr_cluster}}{Matrix of size (\code{nrep}, J), with J = \code{length(unique(group))}. 
#'    Each row is a posterior sample of the distributional cluster allocation variables \eqn{(S_1,\dots,S_J)}. }
#' \item{\code{pi}}{Matrix of size (\code{nrep}, \code{maxK}). 
#'    Each row is a posterior sample of the distributional cluster probabilities \eqn{(\pi_1,\dots,\pi_{maxK})}.}
#' \item{\code{omega}}{3-d array of size (\code{maxL}, \code{maxK}, \code{nrep}).
#'    Each slice is a posterior sample of the observational cluster probabilities. 
#'    In each slice, each column \eqn{k} is a vector (of length \code{maxL}) observational cluster probabilities
#'    \eqn{(\omega_{1,k},\dots,\omega_{L,k})} for distributional cluster \eqn{k}. }
#' \item{\code{alpha}}{Vector of length \code{nrep} of posterior samples of the parameter \eqn{\alpha}.}
#' \item{\code{beta}}{Vector of length \code{nrep} of posterior samples of the parameter \eqn{\beta}.}
#' \item{\code{maxK}}{Vector of length \code{nrep} of the number of distributional DP components used by the slice sampler.}
#' \item{\code{L_plus}}{Vector of length \code{nrep} of posterior samples of the number of observational clusters.}
#' \item{\code{L}}{Vector of length \code{nrep} of posterior samples of the observational Dirichlet dimension.}
#' }
#'
#' @examples 
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' plot(density(y[g==1]), xlim = c(-5,10))
#' lines(density(y[g==2]), col = 2)
#' out <- sample_fiSAN(nrep = 500, y = y, group = g, 
#'                     nclus_start = 2,
#'                     maxK = 20, maxL = 20,
#'                     beta = 1)
#' out 
#' 
#' @references 
#' Frühwirth-Schnatter, S., Malsiner-Walli, G. and Grün, B. (2021).
#' Generalized mixtures of finite mixtures and telescoping sampling. \emph{Bayesian Analysis}, 16(4), 1279–1307. <doi:10.1214/21-BA1294>
#' 
#' Kalli, M., Griffin, J.E., and Walker, S.G. (2011). Slice Sampling Mixture Models, 
#' \emph{Statistics and Computing}, 21, 93–105. <doi:10.1007/s11222-009-9150-y>
#' 
#' Sethuraman, A.J. (1994). A Constructive Definition of Dirichlet Priors, \emph{Statistica Sinica}, 4, 639–650.
#'
#' @export sample_fiSAN
#' @importFrom stats cor var dist hclust cutree rgamma 
sample_fiSAN <- function(nrep, y, group, 
                       maxK = 50, maxL = 50, 
                       m0 = 0, tau0 = 0.1, lambda0 = 3, gamma0 = 2,
                       hyp_alpha1 = 1, hyp_alpha2 = 1,
                       hyp_beta1 = 6, hyp_beta2 = 3, 
                       eps_beta = NULL,
                       alpha = NULL, beta = NULL,
                       warmstart = TRUE, nclus_start = NULL,
                       mu_start = NULL, sigma2_start = NULL,
                       M_start = NULL, S_start = NULL,
                       alpha_start = NULL, beta_start = NULL,
                       progress = TRUE,
                       seed = NULL)
{
  group <- .relabell(group) - 1 
  
  if(is.null(seed)){seed <- round(stats::runif(1,1,10000))}
  set.seed(seed)
  
  params <- list(nrep = nrep, 
                y = y,
                group = group+1, 
                maxK = maxK, 
                maxL = maxL, 
                m0 = m0, tau0 = tau0,
                lambda0 = lambda0, gamma0 = gamma0,
                seed = seed)
  
  if(!is.null(alpha)) { params$alpha <- alpha }
  if(!is.null(beta)) { params$beta <- beta }
  if(is.null(alpha)) { params$hyp_alpha1 <- hyp_alpha1 }
  if(is.null(alpha)) { params$hyp_alpha2 <- hyp_alpha2 }
  if(is.null(beta)) { params$hyp_beta1 <- hyp_beta1 }
  if(is.null(beta)) { params$hyp_beta2 <- hyp_beta2 }
  
  if(is.null(S_start)) { S_start <- rep(0,length(unique(group))) }
  
  # if the initial cluster allocation is passed
  if(!is.null(M_start)) { 
    warmstart <- FALSE
    M_start <- .relabell(M_start)
    
    # and the mean is passed or the variance is passed don't do anything
    
    # if the mean is not passed
    if(is.null(mu_start)) { 
      mu_start <- rep(0,maxL)
      ncl0 <- length(unique(M_start))
      for(l in unique(M_start)) {
        mu_start[l] <- mean(y[M_start == l])
      }
    }
    # if the variance is not passed
    if(is.null(sigma2_start)) { 
      sigma2_start <- rep(0.001,maxL)
      ncl0 <- length(unique(M_start))
      for(l in unique(M_start)) {
        sigma2_start[l] <- var(y[M_start == l])
      }
    }
  } else {
    # if the initial cluster allocation is not passed
    # and you don't want a warmstart
    if(!warmstart){
      M_start <- rep(1, length(y))#sample(0:(maxL-2), length(y), replace = TRUE)
      mu_start <- rep(0, maxL)
      mu_start[1] <- mean(y)
      sigma2_start <- rep(0.001, maxL)
      sigma2_start[1] <- var(y)/2
    }
    
    # if the initial cluster allocation is not passed
    # and you want a warmstart
    if(warmstart){
      mu_start <- rep(0,maxL)
      sigma2_start <- rep(0.001,maxL)
      
      if(is.null(nclus_start)) { nclus_start <- min(c(maxL, 30))}
      M_start <- stats::kmeans(y,
                              centers <- nclus_start, 
                              algorithm="MacQueen",
                              iter.max = 50)$cluster 
      
      nclus_start <- length(unique(M_start))
      mu_start[1:nclus_start] <- sapply(1:nclus_start, function(x) mean(y[M_start == x])) 
      sigma2_start[1:nclus_start] <- sapply(1:nclus_start, function(x) var(y[M_start == x])) 
      sigma2_start[1:nclus_start][sigma2_start[1:nclus_start]==0] <- 0.001
      sigma2_start[is.na(sigma2_start)] <- 0.001
    }
  }
  M_start <- M_start-1
  sigma2_start[is.na(sigma2_start)] <- 0.001
  
  if(is.null(alpha_start)) { alpha_start <- rgamma(1, hyp_alpha1, hyp_alpha2) }
  if(is.null(beta_start)) { beta_start <- rgamma(1, hyp_beta1, hyp_beta2) }
  
  fixed_alpha <- F
  fixed_beta <- F
  if(!is.null(alpha) ) {
    fixed_alpha <- T ; } else { alpha <- 1 }
  if(!is.null(beta) ) {
    fixed_beta <- T ; eps_beta <- 1 } else { beta <- 1 }
  
  if((fixed_beta == F)&(is.null(eps_beta))) {stop("Missing eps parameter for MH step on beta Please provide 'eps_beta' or a fixed 'beta' value.")}
  
  start <- Sys.time()
  out <- sample_ficam_arma(nrep, y, group, 
                           maxK, maxL, 
                           m0, tau0, 
                           lambda0, gamma0,
                           fixed_alpha, fixed_beta,
                           alpha, beta,
                           hyp_alpha1, hyp_alpha2,
                           hyp_beta1, hyp_beta2,
                           mu_start, sigma2_start,
                           M_start, S_start, 
                           alpha_start, beta_start,
                           eps_beta,
                           progress)
  end <- Sys.time()
  
  warnings <- out$warnings
  out[12] <- NULL
  
  out$distr_cluster <- out$distr_cluster + 1 
  out$obs_cluster <- out$obs_cluster + 1 
  
  
  if(length(warnings) == 2) {
    output <- list( "model" = "fiSAN",
                   "params" = params,
                   "sim" = out,
                   "time" = end - start,
                   "warnings" = warnings)
    warning("Increase maxL and maxK: all the provided mixture components were used. Check $warnings to see when it happened.")
  } else if (length(warnings) == 1) {
    if((length(warnings$top_maxK)>0) & (length(warnings$top_maxL)==0)) {
      output <- list( "model" = "fiSAN",
                     "params" = params,
                     "sim" = out,
                     "time" = end - start,
                     "warnings" = warnings)
      warning("Increase maxK: all the provided distributional mixture components were used. Check '$warnings' to see when it happened.")
    }
    
    if((length(warnings$top_maxK)==0) & (length(warnings$top_maxL)>0)) {
      output <- list( "model" = "fiSAN",
                     "params" = params,
                     "sim" = out,
                     "time" = end - start,
                     "warnings" = warnings)
      warning("Increase maxL: all the provided observational mixture components were used. Check '$warnings' to see when it happened.")
    }
  } else {
    output <- list( "model" = "fiSAN",
                   "params" = params,
                   "sim" = out,
                   "time" = end - start )
  }
  
  structure(output, class = c("SANmcmc", class(output)))
  
}
