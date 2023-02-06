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
NLSS_sum = function(res, th=0.95, nstart = 1, nend = 1){

  state0 = res$state0
  states = res$states

  n = dim(res$A)[1]

  n_itrs = nend - nstart + 1

  p = ncol(res$S)
  q = nrow(res$S)

  # tmp = factor(res$S[,,nstart:nend],levels=res$states)
  #
  # dim(tmp) = dim(res$S[,,nstart:nend])

  tmp = res$S[,,nstart:nend]

  beta_coef = apply(tmp,c(1,2),count_rate,states,n_itrs)

  S0 = states[apply(beta_coef,c(2,3),which.max)]

  dim(S0) = dim(res$S[,,1])

  S = S0

  pip = 1 - beta_coef[which(states==state0),,]

  if(th==1){
    th=1-1e-6
  }

  for(i in 1:q){
    for(j in 1:p){
      if(S[i,j]!=(state0) ){
        ind = which(states==S[i,j])
        if(beta_coef[ind,i,j]<=th){
          S[i,j] = state0
        }
      }
    }
  }

  loglik0 = rep(0,n_itrs)

  for(i in 1:n_itrs){
    A1 = res$A[,,nstart+i-1]
    S1 = res$S[,,nstart+i-1]

    loglik0[i] = NLSS_logLik_noise0(res$X,S1,A1,res$K)
  }

  Amean = apply(res$A[,,nstart:nend],c(1,2),mean)

  loglik = NLSS_logLik_noise0(res$X,S,Amean,res$K)
  BIC = -2 * loglik + log(dim(res$X)[1]*dim(res$X)[2])*sum(S!=state0)

  return(list(A=Amean, S_mode =S0, S = S, Z=res$Y, pip = pip,
              BIC=BIC, loglik = loglik, loglik_mcmc = loglik0))

}



count_rate = function(x,states,nitr){
  return( table(factor(x,levels=states))/nitr )
}


