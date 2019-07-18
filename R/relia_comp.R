#' Title
#'
#' @param Slist a list of matrices
#'
#' @return
#' @export
#'
#' @examples
relia_comp = function(Slist){

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

  return(list(r_out = r_out, r_out2 = r_out2) )
}
