#' @title MCMC sampling for NLSS
#' @description The function conducts MCMC sampling for the NLSS model. The number of latent sources need to be specified.
#'
#' @param data a n*p matrix with discrete values, n is the sample size, p is the num of node pairs.
#' @param q the number of latent sources.
#' @param total_iter the number of total iterations.
#' @param burn_in the number of iterations to be discarded as burn-in.
#' @param show_step the frequency for printing the current number of iterations.
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
NLSS = function(data, q=2, kk=1, group_node=NULL,init = list(S=NULL,A=NULL,beta=NULL), total_iter = 1000, burn_in = 500,thin =10, show_step = 100 ){

  n = nrow(data)
  p = ncol(data)
  mindata = min(data)
  up = 0
  while(mindata<=0){
    mindata = mindata + 1
    up = up + 1
  }

  X0 = data + up
  K = max(X0)

  if(is.null(init$A)){
    A0 = matrix(1.0/(q+1),nrow=n,ncol=q+1)
  }
  else{
    A0 = init$A
  }

  if(is.null(init$Y)){
    Y0 = matrix(sample(q+1,n*p,replace = TRUE),nrow = n, ncol = p)
  }
  else{
    Y0 = init$Y0
  }

  if(is.null(init$S)){
    S0 = matrix(0,nrow = q,ncol =p)
  }
  else{
    S0 = init$S
  }

  if(is.null(group_node)){
    group= rep(0,p)
  }
  else{
    L= length(group_node)
    groupmat = matrix(0,L,L)
    tag = 0
    for(i in unique(group_node)){
      for(j in unique(group_node)){
        if(j>=i){
          groupmat[group_node==i,group_node==j] = tag
          groupmat[group_node==j,group_node==i] = tag
          tag = tag + 1
        }
      }
    }
    group = vec_mat(groupmat)
  }

  G = max(group) + 1

  if(is.null(init$beta)){
    beta0 = matrix(1/K, nrow = G*q, ncol= K)
  }else{
    beta0 = init$beta
  }

  gamma0 = 0.5
  alpha0 = c(0.5,0.5)
  res = NLSS_gibbs_sampler_n(X=X0, A0 =A0, S0 = S0, Y0=Y0, beta0 = beta0,group = group,gamma = gamma0,alpha = alpha0 , kk=kk,total_iter = total_iter, burn_in = burn_in, thin =thin, show_step = show_step)

  res$group = group
  #sum_res = NLSS_sum(res)

  #out = list()
  #out$res = res
  #out$sum_res = sum_res

  return(res)
}



