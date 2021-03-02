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
relia_rows = function(Slist, conn_only=TRUE){

  S0 = Slist[[1]]
  if(class(S0)!="matrix"){
    S0 = t(as.matrix(S0))
  }

  r_out = rep(0,nrow(S0))

  for(i in 1:nrow(S0)){
    tmp0 = 0
    if(conn_only){
      ind0 = S0[i,]!=0
    }
    else{
      ind0 = 1:ncol(S0)
    }

    for(j in 2:length(Slist) ){
      S1 = Slist[[j]]
      if(class(S1)!="matrix"){
        S1 = t(as.matrix(S1))
      }
      tmp0 = tmp0 + mean(S0[i,ind0]==S1[i,ind0])
    }

    tmp1 = tmp0

    for(j in 2:length(Slist) ){
      for(k in 1:nrow(S0)){
        if(k!=i) {
          S1 = Slist[[j]]
          if(class(S1)=="numeric"){
            S1 = t(as.matrix(S1))
          }
          tmp1 = tmp1 + mean(S0[i,ind0]==S1[k,ind0])
        }
      }
    }
    tmp0 = tmp0 / (length(Slist)-1)
    tmp1 = tmp1/  (length(Slist)-1) / nrow(S0)

    r_out[i] = (tmp0 - tmp1) / (1 - tmp1)
  }

  return(r = r_out )
}

#' @export
relia_rows_bygroup = function(Slist,group,conn_only=TRUE){

  S0 = Slist[[1]]
  if(class(S0)!="matrix"){
    S0 = t(as.matrix(S0))
  }

  L = length(group)

  groupmat = matrix(0,L,L)
  tag = 1
  for(i in unique(group)){
    for(j in unique(group)){
      if(j>=i){
        groupmat[group==i,group==j] = tag
        groupmat[group==j,group==i] = tag
        tag = tag + 1
      }
    }
  }
  group_vec = vec_mat(groupmat)

  ind = unique(group_vec)
  L = length(ind)

  r_out = matrix(0,nrow=L, ncol=nrow(S0))
  r_out2 = matrix(0,nrow=L, ncol=nrow(S0))

  tag = 1
  for(g in ind){

    S0_g = S0[,group_vec==g]

    for(i in 1:nrow(S0)){
      tmp0 = 0

      if(conn_only){
        ind0 = (S0_g[i,]!=0)
      }
      else{
        ind0 = 1:ncol(S0_g)
      }


      for(j in 2:length(Slist) ){
        S1 = Slist[[j]]
        if(class(S1)!="matrix"){
          S1 = t(as.matrix(S1))
        }
        S1_g = S1[,group_vec==g]
        tmp0 = tmp0 + mean(S0_g[i,ind0]==S1_g[i,ind0])
      }

      tmp1 = tmp0

      for(j in 2:length(Slist) ){
        for(k in 1:nrow(S0)){
          if(k!=i) {
            S1 = Slist[[j]]
            if(class(S1)=="numeric"){
              S1 = t(as.matrix(S1))
            }
            S1_g = S1[,group_vec==g]
            tmp1 = tmp1 + mean(S0_g[i,ind0]==S1_g[k,ind0])
          }
        }
      }
      tmp0 = tmp0 / (length(Slist)-1)
      tmp1 = tmp1/  (length(Slist)-1) / nrow(S0)

      r_out[g,i] = (tmp0 - tmp1) / (1 - tmp1)
      r_out2[g,i] = tmp0
    }


  }


  out = list()
  out$r = r_out
  out$group = group_vec
  return( out)
}

