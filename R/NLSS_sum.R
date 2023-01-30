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
NLSS_sum = function(res, th=0.95, k0 = 0, nstart = 1, nend = 1){

  n = dim(res$A)[1]

  n_itrs = nend - nstart + 1

  p = ncol(res$S)
  q = nrow(res$S)

  tmp = factor(res$S[,,nstart:nend],levels=1:res$K)

  dim(tmp) = dim(res$S[,,nstart:nend])

  beta_coef = apply(tmp,c(1,2),count_rate,res$K,n_itrs)

  S0 = apply(beta_coef,c(2,3),which.max)

  S = S0

  pip = 1 - beta_coef[k0,,]

  if(th==1){
    th=1-1e-6
  }

  for(i in 1:q){
    for(j in 1:p){
      if(S[i,j]!=(k0) ){
        if(beta_coef[S[i,j],i,j]<=th){
          S[i,j] = k0
        }
      }
    }
  }

  loglik0 = rep(0,n_itrs)

  for(i in 1:n_itrs){
    A1 = res$A[,,nstart+i-1]
    S1 = res$S[,,nstart+i-1]

    loglik0[i] = DLSS_logLik_noise0(res$X,S1,A1,res$K)
  }

  Amean = apply(res$A[,,nstart:nend],c(1,2),mean)

  loglik = DLSS_logLik_noise0(res$X,S,Amean,res$K)
  BIC = -2 * loglik_S + log(dim(res$X)[1]*dim(res$X)[2])*sum(S!=k0)

  return(list(A=Amean, S_mode =S0, S = S, Z=res$Y, pip = pip,
              BIC=BIC, loglik = loglik, loglik_mcmc = loglik0))

}


count_rate = function(x,K,nitr){
  return( table(factor(x,levels=1:K))/nitr )
}
