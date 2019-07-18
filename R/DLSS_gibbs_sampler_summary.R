#' Title
#'
#' @param res result from the function NLSS
#'
#' @return
#' @export
#'
#' @examples
DLSS_gibbs_sampler_summary = function(res){

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
