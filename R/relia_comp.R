#' @title Reliability of sources
#' @description The function calculates the reliability index of sources for NLSS. The input is a list of
#' matrices with the first matrix being estimated sources from the whole dataset and the following matrices
#' being estimated sources from bootsrap samples of the whole dataset.
#'
#' @param Slist a list of matrices with rows representing vectorization of latent sources.
#'
#' @return Reliability index.
#' @export
#'
#' @examples
relia_rows = function(Slist){

  S0 = Slist[[1]]
  if(class(S0)!="matrix"){
    S0 = t(as.matrix(S0))
  }

  r_out = rep(0,nrow(S0))
  r_out2 = rep(0,nrow(S0))

  for(i in 1:nrow(S0)){
    tmp0 = 0
    for(j in 2:length(Slist) ){
      S1 = Slist[[j]]
      if(class(S1)!="matrix"){
        S1 = t(as.matrix(S1))
      }
      tmp0 = tmp0 + mean(S0[i,]==S1[i,])
    }

    tmp1 = tmp0

    for(j in 2:length(Slist) ){
      for(k in 1:nrow(S0)){
        if(k!=i) {
          S1 = Slist[[j]]
          if(class(S1)=="numeric"){
            S1 = t(as.matrix(S1))
          }
          tmp1 = tmp1 + mean(S0[i,]==S1[k,])
        }
      }
    }
    tmp0 = tmp0 / (length(Slist)-1)
    tmp1 = tmp1/  (length(Slist)-1) / nrow(S0)

    r_out[i] = (tmp0 - tmp1) / (1 - tmp1)
    r_out2[i] = tmp0
  }

  return(r = r_out )
}
