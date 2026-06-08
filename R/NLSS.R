#' @title MCMC Sampling for the NLSS Model
#' This function performs Markov chain Monte Carlo (MCMC) sampling for the
#' Network Latent Source Separation (NLSS) model. Given a discrete-valued
#' network data matrix, the function estimates the latent source matrix, the
#' mixing coefficient matrix, and related posterior samples.
#'
#' @param data An \eqn{n \times p} matrix of discrete-valued connection states,
#'   where \eqn{n} is the sample size and \eqn{p} is the number of node pairs.
#' @param states A vector containing all possible discrete states in the observed
#'   data and latent sources. The default is \code{min(data):max(data)}.
#' @param state0 The background or null state. The default is \code{min(data)}.
#' @param q An integer specifying the number of latent sources.
#' @param q0 An integer specifying the maximum number of non-background sources
#'   allowed for each node pair. The default is \code{q}.
#' @param init A list of initial values for the MCMC sampler. It may contain
#'   \code{S}, \code{A}, and \code{Y}, corresponding to the initial latent source
#'   matrix, mixing coefficient matrix, and latent allocation matrix,
#'   respectively. Missing components are initialized automatically.
#' @param total_iter An integer specifying the total number of MCMC iterations.
#'   The default is \code{1000}.
#' @param burn_in An integer specifying the number of initial MCMC iterations to
#'   be discarded as burn-in. The default is \code{500}.
#' @param thin An integer specifying the thinning interval for saving posterior
#'   samples. The default is \code{10}.
#' @param show_step An integer specifying how often the current iteration number
#'   is printed during MCMC sampling. The default is \code{100}.
#' @param joint A logical value indicating whether joint updating is used in the
#'   MCMC sampler. The default is \code{TRUE}.
#'
#' @return A list containing posterior samples and model information from the
#'   MCMC sampler. The returned object includes, but is not necessarily limited
#'   to, the following components:
#' \describe{
#'   \item{\code{A}}{Posterior samples of the mixing coefficient matrix.}
#'   \item{\code{S}}{Posterior samples of the latent source matrix.}
#'   \item{\code{Y}}{Posterior samples of the latent allocation matrix.}
#'   \item{\code{X}}{The input discrete-valued data matrix.}
#'   \item{\code{K}}{The number of possible discrete states.}
#'   \item{\code{states}}{The vector of possible discrete states.}
#'   \item{\code{state0}}{The background or null state.}
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
NLSS = function(data,  q=2, q0 = q, init = list(S=NULL,A=NULL,Y=NULL),
                states = (min(data):max(data)), state0 = min(data),
                total_iter = 1000, burn_in = 500, thin = 10, show_step = 100,
                joint=TRUE){

  n = nrow(data)
  p = ncol(data)

  K = length(states)

  Sspace = permutations(K,q,v=states, repeats.allowed=TRUE)
  n_nonzero = apply( (Sspace!=(state0) ),1,sum)
  Sspace = Sspace[ which(n_nonzero<=q0), ]

  stat_link = matrix(0,nrow = nrow(Sspace), ncol = K)

  index = matrix(0,nrow = 1, ncol =q)

  res_index = sum_subset2(index,q,1,q0)

  for(k in 1:K){
    s0 = states[k]
    tmp = (Sspace[,1:q]==s0)
    for(i in 1:nrow(stat_link)){
      comb = NULL
      if(s0!=state0){
        for(j in 1:q ){
          if(tmp[i,j]==1){
            comb = paste0(comb,"a",j)
          }
        }
        if(!is.null(comb)){
          stat_link[i,k] = which(res_index$eqns==comb)
        }
        else{
          stat_link[i,k] = 1
        }
      }
      else{
        for(j in 1:q ){
          if(tmp[i,j]==0){
            comb = paste0(comb,"a",j)
          }
        }
        if(!is.null(comb)){
          stat_link[i,k] = which(res_index$eqns==comb) + length(res_index$eqns)
        }
        else{
          stat_link[i,k] = 1 + length(res_index$eqns)
        }
      }
    }
  }

  if(is.null(init$A)){
    A0 = matrix(1.0/(q+1),nrow=n,ncol=q+1)
  }
  else{
    A0 = init$A
  }

  if(is.null(init$Y)){
    Y0 = matrix(sample(q,n*p,replace = TRUE),nrow = n, ncol = p)
  }
  else{
    Y0 = init$Y
  }

  if(is.null(init$S)){
    S0 = matrix(state0,nrow = q,ncol =p)
  }
  else{
    S0 = init$S
  }

  alpha0 = rep(0.5,q+K)

  res = NLSS_gibbs_sampler(X=data, A0 =A0, S0 = S0, Y0=Y0,
                           stat_link = stat_link,
                           joint=joint, alpha = alpha0, Sspace=Sspace,
                           states=states,
                           q0 = q0, K=K,
                           total_iter = total_iter, burn_in = burn_in,
                           thin = thin, show_step = show_step)

  res$states = states
  res$state0 = state0

  return(res)
}

