#' @title Summary of the MCMC result for NLSS
#' @description The function summarizes the MCMC result and returns the posterior mean
#' of A, the posterior mode of S, beta coefficient (frequency of each discrete value of S
#' among the MCMC samples) and the log-likelihood trace.
#'
#' @param res result from the function NLSS
#'
#' @return
#'
#' @examples
NLSS_sum = function(res){

  n = dim(res$A)[1]
  save_iters= dim(res$S)[3]
  p = ncol(res$S)
  q = nrow(res$S)
  temp = factor(res$S,levels=1:res$K)
  dim(temp) = dim(res$S)
  beta_coef = array(0,dim=c(q,p,res$K))
  for(j in 1:p){
    for(l in 1:q)
      beta_coef[l,j,] = table(factor(temp[l,j,],levels=1:res$K))/save_iters;
  }

  A_coef = apply(res$A,c(1,2),mean)

  S = apply(beta_coef,c(1,2),which.max)

  loglik0 = rep(0,save_iters)

  for(i in 1:save_iters){
    A1 = res$A[,,i]
    loglik0[i] = DLSS_logLik_noise(res$X,S,A1,res$K)
  }

  return(list(A=A_coef,beta_coef = beta_coef,S = S, logLik = loglik0))

}
