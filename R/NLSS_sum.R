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
#' @export
NLSS_sum = function(res){

  mindata = min(res$X)
  up = 0
  while(mindata<=0){
    mindata = mindata + 1
    up = up + 1
  }

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
  S = S - up

  beta = apply(res$beta,c(1,2),mean)

  loglik0 = rep(0,save_iters)
  G0 = max(res$group)+1

  for(i in 1:save_iters){
    A1 = res$A[,,i]
    beta1 = res$beta[,,i]
    loglik0[i] = NLSS_logLik_noise(res$X,A1,beta1,res$group,res$K,G0)
  }

  return(list(A=A_coef,beta=beta,beta_coef = beta_coef, S = S, Y=res$Y, logLik = loglik0 ))

}

#' Title
#'
#' @param X
#' @param res
#' @param res_sum
#' @param start
#' @param end
#' @param K
#'
#' @return
#' @export
#'
#' @examples
log_lik_A = function(X,res,res_sum,start=1,end=2,K){

  out = rep(0,end-start+1)
  tag = 1
  #A1 = res_sum$A
  G0 = max(res$group)+1
  for(i in start:end){
    S1 = res$S[,,i]
    beta1 = res$beta[,,i]
    A1 = res$A[,,i]
    out[tag] = DLSS_logLik_noise0(X,S1,A1,beta1,K)
    tag = tag + 1
  }
  return(out)
}

#' @export
log_lik_group = function(X,res,res_sum,g0,start=1,end=2,K){

  out = rep(0,end-start+1)
  tag = 1
  A1 = res_sum$A
  G0 = max(res$group)+1
  for(i in start:end){
    S1 = res$S[,,i]
    beta1 = res$beta[,,i]
    #A1 = res$A[,,i]
    #out[tag] = NLSS_logLik_noise_group(X,A1,beta1,res$group,g0,K,G0)
    out[tag] = DLSS_logLik_noise0_group(X,S1,A1,res$group,g0,K,G0)
    tag = tag + 1
  }
  return(out)
}



#' @export
log_lik_A_beta = function(X,res,res_sum,start=1,end=2,K){

  out = rep(0,end-start+1)
  tag = 1
  #A1 = res_sum$A
  G0 = max(res$group)+1
  for(i in start:end){
    S1 = res$S[,,i]
    beta1 = res$beta[,,i]
    A1 = res$A[,,i]
    out[tag] = NLSS_logLik_noise(X,A1,beta1,res$group,K,G0)
    tag = tag + 1
  }
  return(out)
}
