#' @title MCMC sampling for NLSS
#' @description The function conducts MCMC sampling for the NLSS model.
#'
#' @param data a n*p matrix with discrete-valued connection states, n is the sample size, p is the num of node pairs.
#' @param k0 the zero states.
#' @param q the number of latent sources.
#' @param init a list of initial values for S and A.
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
NLSS = function(data, k0 = 1, q=2,  init = list(S=NULL,A=NULL,beta=NULL),
                total_iter = 1000, burn_in = 500,thin =10, show_step = 100, sprs = 2,
                joint=TRUE ){

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

  Sspace = permutations(K,q,repeats.allowed=TRUE)
  n_nonzero = apply( (Sspace!=(k0+up ) ),1,sum)
  Sspace = Sspace[ which(n_nonzero<=sprs), ]

  stat_link = matrix(0,nrow = nrow(Sspace), ncol = K)

  index = matrix(0,nrow = 1, ncol =q)

  res_index = sum_subset2(index,q,1,sprs)
  #res_index$eqns = c(res_index$eqns,res_index$eqns)

  for(k in 1:K){
    tmp = (Sspace[,1:q]==k)
    for(i in 1:nrow(stat_link)){
      comb = NULL
      if(k!=(k0+up)){
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

  if(is.null(init$B)){
    B0 = matrix(1.0/(q+1)/K,nrow=n,ncol=K)
  }
  else{
    B0 = init$B
  }

  if(is.null(init$Y)){
    Y0 = matrix(sample(q+K,n*p,replace = TRUE),nrow = n, ncol = p)
  }
  else{
    Y0 = init$Y0
  }

  if(is.null(init$S)){
    S0 = matrix(2,nrow = q,ncol =p)
  }
  else{
    S0 = init$S
  }

  alpha0 = rep(0.5,q+K)

  res = NLSS_gibbs_sampler(X=X0, A0 =A0, S0 = S0, Y0=Y0, stat_link = stat_link,
                           joint=joint, alpha = alpha0, Sspace=Sspace,
                           sprs = sprs, K=K,
                           total_iter = total_iter, burn_in = burn_in, thin = thin, show_step = show_step)

  res$stat_link = stat_link
  res$Sspace = Sspace
  res$res_index = res_index

  return(res)
}


