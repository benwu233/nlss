#' @title MCMC sampling for NLSS
#' @description The function conducts MCMC sampling for a NLSS model. The number of latent sources need to be specified.
#'
#' @param data a n*p matrix with discrete value entries, n is the sample size, p is the num of node pairs,
#' @param q the number of latent sources.
#' @param total_iter the number of total iterations.
#' @param burn_in the number of iterations to be discarded as burn-in.
#' @param show_step the number of iterations as a gap between two consecutive prints for monitoring the algorithm.
#'
#' @return a R list containing the following terms:
#' \describe{
#'   \item{A}{the estimated mixing coefficent matrix.}
#'   \item{beta_coef}{the frequency of each discrete value of S appearing
#' in the MCMC samples.}
#'   \item{S}{the estimated latent source matrix.}
#'   \item{logLik}{the log-likelihood trace.}
#' }
#'
#' @importFrom RcppParallel RcppParallelLibs
#' @import Rcpp
#' @import MCMCpack
#' @useDynLib nlss
#'
#' @export
#'
#' @examples
NLSS = function(data, q=2, total_iter = 1000, burn_in = 500, show_step = 100 ){

  q = q + 1
  n = nrow(data)
  p = ncol(data)
  A0=matrix(1.0/q,nrow=n,ncol=q)
  mindata = min(data)
  up = 0

  while(mindata<=0){
    mindata = mindata + 1
    up = up + 1
  }

  X0 = data + up
  gamma0 = as.numeric(table(X0))
  gamma0 = gamma0*10/sum(gamma0)
  alpha0 = c(0.5,0.1)
  res = DLSS_gibbs_sampler_n(X=X0, A0 =A0,gamma = gamma0,alpha = alpha0 ,total_iter = total_iter, burn_in = burn_in, show_step = show_step)
  sum_res = NLSS_sum(res)

  sum_res$S = sum_res$S - up

  return(sum_res)
}
