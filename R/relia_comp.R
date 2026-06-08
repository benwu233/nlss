#' @title Reliability of Latent Sources
#' @description
#' This function calculates the reliability index of latent sources estimated by
#' the NLSS model. The input is a list of source matrices, where the first matrix
#' is the source estimate obtained from the whole dataset and the remaining
#' matrices are source estimates obtained from bootstrap samples. The reliability
#' index compares the reproducibility of the same source across bootstrap samples
#' with the reproducibility expected from other sources.
#'
#' @param Slist A list of source matrices. Each matrix should have rows
#'   representing vectorized latent sources and columns representing node pairs
#'   or features. The first matrix is treated as the reference source estimate,
#'   while the remaining matrices are treated as bootstrap estimates. All matrices
#'   should have the same dimensions and their rows should already be matched.
#' @param conn_only A logical value indicating whether reliability should be
#'   computed only over non-background connections in the reference source. If
#'   \code{TRUE}, only entries of the reference source that are not equal to
#'   zero are used. If \code{FALSE}, all entries are used. The default is
#'   \code{TRUE}.
#' @return A numeric vector containing the reliability index for each latent
#'   source. Each element corresponds to one row of the reference source matrix.
#' @export
#'
relia_rows = function(Slist, conn_only=TRUE){

  S0 = Slist[[1]]
  if(class(S0)[1]!="matrix"){
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
      if(class(S1)[1]!="matrix"){
        S1 = t(as.matrix(S1))
      }
      tmp0 = tmp0 + mean(S0[i,ind0]==S1[i,ind0])
    }

    tmp1 = tmp0

    for(j in 2:length(Slist) ){
      for(k in 1:nrow(S0)){
        if(k!=i) {
          S1 = Slist[[j]]
          if(class(S1)[1]=="numeric"){
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


relia_idx = function(Slist, group, g0,g1){

  S0 = Slist[[1]]
  if(class(S0)[1]!="matrix"){
    S0 = t(as.matrix(S0))
  }

  S2 = matrix(0,236,236)
  S2[group==g0,group==g1] = 1
  S2[group==g1,group==g0] = 1
  S3 = vec_mat(S2)


  r_out = rep(0,nrow(S0))

  for(i in 1:nrow(S0)){
    tmp0 = 0
    ind0 = (S3 ==1)

    for(j in 2:length(Slist) ){
      S1 = Slist[[j]]
      if(class(S1)[1]!="matrix"){
        S1 = t(as.matrix(S1))
      }
      tmp0 = tmp0 + mean(S0[i,ind0]==S1[i,ind0])
    }

    tmp1 = tmp0

    for(j in 2:length(Slist) ){
      for(k in 1:nrow(S0)){
        if(k!=i) {
          S1 = Slist[[j]]
          if(class(S1)[1]=="numeric"){
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

#' @title Reliability of Latent Sources by Node Group
#'
#' @description
#' This function calculates the reliability index of latent sources separately
#' for each within-group or between-group block of network connections. The first
#' matrix in \code{Slist} is treated as the reference estimate from the whole
#' dataset, and the remaining matrices are treated as bootstrap estimates.
#'
#' @param Slist A list of source matrices. Each matrix should have rows
#'   representing vectorized latent sources and columns representing node pairs.
#'   The first matrix is treated as the reference source estimate, while the
#'   remaining matrices are treated as bootstrap estimates. All matrices should
#'   have the same dimensions and their rows should already be matched.
#' @param group A vector of group labels for network nodes. Its length should be
#'   equal to the number of nodes in the network.
#' @param conn_only A logical value indicating whether reliability should be
#'   computed only over non-background connections in the reference source within
#'   each group block. If \code{TRUE}, only entries of the reference source that
#'   are not equal to zero are used. If \code{FALSE}, all entries are used. The
#'   default is \code{TRUE}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{r}}{A matrix of reliability indices. Rows correspond to
#'   within-group or between-group connection blocks, and columns correspond to
#'   latent sources.}
#'   \item{\code{group}}{A vector indicating the group-block membership of each
#'   vectorized node pair. This vector has the same length as the number of
#'   columns in the source matrices.}
#' }
#'
#' @export
relia_rows_bygroup = function(Slist,group,conn_only=TRUE){

  S0 = Slist[[1]]
  if(class(S0)[1]!="matrix"){
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
        if(class(S1)[1]!="matrix"){
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
            if(class(S1)[1]=="numeric"){
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

