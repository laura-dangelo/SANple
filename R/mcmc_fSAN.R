#' Sample fSAN with sparse mixtures
#' 
#' @description \code{sample_fSAN} is used to perform posterior inference under the finite shared atoms nested (fSAN) model with Gaussian likelihood (originally proposed in D'Angelo et al., 2023). 
#' The model uses overfitted (sparse) Dirichlet mixtures (Malsiner-Walli et al., 2016) at both the observational and distributional level.
#' 
#' 
#' @usage 
#' sample_fSAN(y, group, ...)
#' 
#' 
#' @param nrep Number of MCMC iterations.
#' @param burn Number of discarded iterations.
#' @param y Vector of observations.
#' @param group Vector of the same length of y indicating the group membership (numeric).
#' @param maxK Maximum number of distributional clusters \eqn{K} (default = 50).
#' @param maxL Maximum number of observational clusters \eqn{L} (default = 50).
#' @param m0,tau0,lambda0,gamma0 Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}. Default is (0, 0.1, 3, 2).
#' @param alpha Distributional Dirichlet parameter if fixed (optional). The distribution is Dirichlet( \code{rep(alpha, maxK)} ).
#' @param beta Observational Dirichlet parameter if fixed (optional). The distribution is Dirichlet( \code{rep(beta, maxL)} ).
#' @param warmstart,nclus_start Initialization of the observational clustering. 
#' \code{warmstart} is logical parameter (default = \code{TRUE}) of whether a kmeans clustering should be used to initialize the chains.
#' An initial guess of the number of observational clusters can be passed via the \code{nclus_start} parameter (optional)
#' @param mu_start,sigma2_start,M_start,S_start Starting points of the MCMC chains (optional). Default is \code{nclus_start = min(c(maxL, 30))}.
#' \code{mu_start, sigma2_start} are vectors of length \code{maxL}. 
#' \code{M_start} is a vector of observational cluster allocation of length N.
#' \code{S_start} is a vector of observational cluster allocation of length J.
#' @param progress show a progress bar? (logical, default TRUE).
#' @param seed set a fixed seed.
#'
#'
#' @details 
#' \strong{Data structure}
#' 
#' The overfitted mixture common atoms model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups. 
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
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,\dots,K\}}, with 
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,\dots,K.}
#' The distribution of the probabilities is \eqn{(\pi_1,\dots,\pi_{K})\sim Dirichlet_K(\alpha,\dots,\alpha)}. 
#' 
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,\dots,K \, ; \: l = 1,\dots,L. }
#' The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim Dirichlet_L(\beta,\dots,\beta)} for all \eqn{k = 1,\dots,K}. 
#'
#'
#'
#' @return \code{sample_fSAN} returns four objects:
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
#'    Each row is a posterior sample of the distributional cluster probabilities \eqn{(\pi_1,\dots,\pi_{K})}.}
#' \item{\code{omega}}{3-d array of size (\code{maxL}, \code{maxK}, \code{nrep}).
#'    Each slice is a posterior sample of the observational cluster probabilities. 
#'    In each slice, each column \eqn{k} is a vector (of length \code{maxL}) observational cluster probabilities
#'    \eqn{(\omega_{1,k},\dots,\omega_{L,k})} for distributional cluster \eqn{k}. }
#' \item{\code{alpha}}{Vector of length \code{nrep} of posterior samples of the parameter \eqn{\alpha}.}
#' \item{\code{beta}}{Vector of length \code{nrep} of posterior samples of the parameter \eqn{\beta}.}
#' }
#'
#' @examples 
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' plot(density(y[g==1]), xlim = c(-5,10))
#' lines(density(y[g==2]), col = 2)
#' out <- sample_fSAN(y = y, group = g)
#' out 
#' 
#' @references D’Angelo, L., Canale, A., Yu, Z., and Guindani, M. (2023). 
#' Bayesian nonparametric analysis for the detection of spikes in noisy calcium imaging data. \emph{Biometrics}, 79(2), 1370–1382. <doi:10.1111/biom.13626>
#' 
#' Malsiner-Walli, G., Frühwirth-Schnatter, S. and Grün, B. (2016). 
#' Model-based clustering based on sparse finite Gaussian mixtures. Statistics and Computing 26, 303–324. <doi:10.1007/s11222-014-9500-2>
#'
#' @export sample_fSAN
#' @importFrom stats cor var dist hclust cutree rgamma 
sample_fSAN <- function(y, group, 
                        prior_param = list(
                                           m0 = 0, 
                                           tau0 = 0.1, 
                                           lambda0 = 3, 
                                           gamma0 = 2,
                                           alpha = NULL, beta = NULL),
                        model_param = list(
                                           nrep = 1000, burn = 500, 
                                           maxK = 50, maxL = 50, 
                                           warmstart = TRUE, nclus_start = NULL,
                                           mu_start = NULL, sigma2_start = NULL, 
                                           M_start = NULL, S_start = NULL,
                                           progress = TRUE,
                                           seed = NULL) )
{
  group <- .relabel(group) - 1 
  
  if(is.null(model_param$seed)){model_param$seed <- round(stats::runif(1,1,10000))}
  set.seed(model_param$seed)
  
  #----------------------------------------------------
  warmstart = model_param$warmstart
  nclus_start = model_param$nclus_start 
  mu_start = model_param$mu_start
  sigma2_start = model_param$sigma2_start
  M_start = model_param$M_start
  S_start = model_param$S_start
  progress = model_param$progress
  seed = model_param$seed         
  #----------------------------------------------------
  
  params <- list(nrep = model_param$nrep, 
                 burn = model_param$burn,
                 y = y,
                 group = group+1, 
                 Nj = tapply(y,group, length),
                 maxK = model_param$maxK, 
                 maxL = model_param$maxL, 
                 m0 =      prior_param$m0, 
                 tau0 =    prior_param$tau0,
                 lambda0 = prior_param$lambda0, 
                 gamma0 =  prior_param$gamma0,
                 seed = seed)
  
  if(!is.null(prior_param$alpha)) { params$alpha <- prior_param$alpha }
  if(!is.null(prior_param$beta)) { params$beta <- prior_param$beta }

  if(is.null(S_start)) { S_start <- rep(0,length(unique(group))) }
  
  # if the initial cluster allocation is passed
  if(!is.null(M_start)) { 
    warmstart <- FALSE
    M_start <- .relabel(M_start)
    
    # and the mean is passed or the variance is passed don't do anything
    
    # if the mean is not passed
    if(is.null(mu_start)) { 
      mu_start <- rep(0,model_param$maxL)
      ncl0 <- length(unique(M_start))
      for(l in unique(M_start)) {
        mu_start[l] <- mean(y[M_start == l])
      }
    }
    # if the variance is not passed
    if(is.null(sigma2_start)) { 
      sigma2_start <- rep(0.001,model_param$maxL)
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
      mu_start <- rep(0, model_param$maxL)
      mu_start[1] <- mean(y)
      sigma2_start <- rep(0.001, model_param$maxL)
      sigma2_start[1] <- var(y)/2
    }
    
    # if the initial cluster allocation is not passed
    # and you want a warmstart
    if(warmstart){
      mu_start <- rep(0,model_param$maxL)
      sigma2_start <- rep(0.001,model_param$maxL)
      
      if(is.null(nclus_start)) { nclus_start <- min(c(model_param$maxL, 30))}
      M_start <- stats::kmeans(y,
                               centers = nclus_start, 
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
  
  if(!is.null(prior_param$alpha) ) {
    fixed_alpha <- T ; 
    alpha_start <- prior_param$alpha
    eps_alpha <- 1
  } else { prior_param$alpha <- 1 }

  if(!is.null(prior_param$beta) ) {
    beta_start <- prior_param$beta
    fixed_beta <- T ; 
    eps_beta <- 1
  } else { prior_param$beta <- 1 }
  
  start <- Sys.time()
  out <- sample_fSAN_cpp(nrep = model_param$nrep, 
                         burn = model_param$burn, 
                         y = y, group = group, 
                         maxK = model_param$maxK, 
                         maxL = model_param$maxL, 
                         m0 = prior_param$m0, 
                         tau0 = prior_param$tau0, 
                         lambda0 = prior_param$lambda0, 
                         gamma0 = prior_param$gamma0,
                         alpha = prior_param$alpha, 
                         beta =  prior_param$beta,
                         mu_start = mu_start, 
                         sigma2_start = sigma2_start,
                         M_start = M_start, 
                         S_start =  S_start, 
                         progressbar = progress)
  end <- Sys.time()
  
  out$distr_cluster <- (out$distr_cluster + 1)
  out$obs_cluster <- (out$obs_cluster + 1)
  
  output <- list( "model" = "fSAN",
                  "params" = params,
                  "sim" = out,
                  "time" = end - start)
  
  structure(output, class = c("SANmcmc",class(output)))
  
}
