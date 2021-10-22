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
NLSS_sum = function(res, nstart = 1, nend = 1){

  mindata = min(res$X)
  up = 0
  while(mindata<=0){
    mindata = mindata + 1
    up = up + 1
  }

  n = dim(res$A)[1]
  #save_iters= dim(res$S)[3]

  n_itrs = nend - nstart + 1

  p = ncol(res$S)
  q = nrow(res$S)

  temp = factor(res$S[,,nstart:nend],levels=1:res$K)

  dim(temp) = dim(res$S[,,nstart:nend])


  # beta_coef = array(0,dim=c(q,p,res$K))
  # for(j in 1:p){
  #   for(l in 1:q)
  #     beta_coef[l,j,] = table(factor(temp[l,j,],levels=1:res$K))/n_itrs;
  # }

  beta_coef = apply(temp,c(1,2),count_rate,res$K,n_itrs)

  A_coef = apply(res$A[,,nstart:nend],c(1,2),mean)

  S = apply(beta_coef,c(2,3),which.max)
  S = S - up

  beta = apply(res$beta[,,nstart:nend],c(1,2),mean)

  loglik0 = rep(0,n_itrs)
  G0 = max(res$group)+1

  for(i in 1:n_itrs){
    A1 = res$A[,,nstart+i-1]
    beta1 = res$beta[,,nstart+i-1]
    loglik0[i] = NLSS_logLik_noise(res$X,A1,beta1,res$group,res$K,G0)
  }

  return(list(A=A_coef,beta=beta,beta_coef = beta_coef, S = S, Y=res$Y, logLik = loglik0 ))

}



#' Title
#'
#' @param res
#' @param nstart
#' @param nend
#'
#' @return
#' @export
#'
#' @examples
NLSS_sum_old = function(res, nstart = 1, nend = 1){

  mindata = min(res$X)
  up = 0
  while(mindata<=0){
    mindata = mindata + 1
    up = up + 1
  }

  n = dim(res$A)[1]
  #save_iters= dim(res$S)[3]

  n_itrs = nend - nstart + 1

  p = ncol(res$S)
  q = nrow(res$S)

  temp = factor(res$S[,,nstart:nend],levels=1:res$K)

  dim(temp) = dim(res$S[,,nstart:nend])

  beta_coef = array(0,dim=c(q,p,res$K))
  for(j in 1:p){
    for(l in 1:q)
      beta_coef[l,j,] = table(factor(temp[l,j,],levels=1:res$K))/n_itrs;
  }

  #beta_coef = apply(temp,c(1,2),count_rate,res$K,n_itrs)

  A_coef = apply(res$A[,,nstart:nend],c(1,2),mean)

  S = apply(beta_coef,c(1,2),which.max)
  S = S - up

  beta = apply(res$beta[,,nstart:nend],c(1,2),mean)

  loglik0 = rep(0,n_itrs)
  G0 = max(res$group)+1

  for(i in 1:n_itrs){
    A1 = res$A[,,nstart+i-1]
    beta1 = res$beta[,,nstart+i-1]
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


count_rate = function(x,K,nitr){
  return( table(factor(x,levels=1:K))/nitr )
}
