#' @title Summary of the MCMC result for NLSS
#' @description The function summarizes the MCMC result returned by \code{NLSS}. It computes
#' the posterior mean of the mixing/loading matrix \code{A}, the element-wise
#' posterior mode of the latent source matrix \code{S}, the posterior inclusion
#' probability of non-background states, the BIC value, and the log-likelihood
#' trace over selected MCMC iterations.
#'
#' @param res A list returned by the function \code{NLSS}. It should contain
#'   MCMC samples of the mixing/loading matrix \code{A} and the latent source
#'   matrix \code{S}, together with other components such as \code{X}, \code{Y},
#'   \code{K}, \code{states}, and \code{state0}.
#' @param th A numeric value between 0 and 1 specifying the posterior frequency
#'   threshold for retaining a non-background state in the estimated latent source
#'   matrix.
#' @param nstart An integer specifying the first MCMC iteration used for
#'   posterior summarization.
#' @param nend An integer specifying the last MCMC iteration used for posterior
#'   summarization.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{A}}{Posterior mean of the mixing/loading matrix \code{A}
#'   over iterations \code{nstart:nend}.}
#'   \item{\code{S_mode}}{Element-wise posterior mode of the latent source
#'   matrix \code{S} over iterations \code{nstart:nend}.}
#'   \item{\code{S}}{Thresholded estimate of the latent source matrix. It is
#'   obtained from \code{S_mode} by replacing weakly supported non-background
#'   states with the background state \code{state0}.}
#'   \item{\code{Z}}{The transformed or reduced data matrix stored as
#'   \code{res$Y}.}
#'   \item{\code{pip}}{Posterior inclusion probability for each entry of
#'   \code{S}, defined as one minus the posterior frequency of the background
#'   state \code{state0}.}
#'   \item{\code{BIC}}{Bayesian information criterion computed from the
#'   thresholded latent source matrix \code{S} and the posterior mean estimate
#'   of \code{A}.}
#'   \item{\code{loglik}}{Log-likelihood evaluated at the posterior mean
#'   estimate of \code{A} and the thresholded estimate of \code{S}.}
#'   \item{\code{loglik_mcmc}}{Log-likelihood trace evaluated at each MCMC
#'   sample from iterations \code{nstart} to \code{nend}.}
#' }
#'
#' @export
NLSS_sum = function(res, th=0.95, nstart = 1, nend = 1){

  state0 = res$state0
  states = res$states

  n = dim(res$A)[1]

  n_itrs = nend - nstart + 1

  p = ncol(res$S)
  q = nrow(res$S)

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


