#' @title Match Rows of Discrete-Valued Source Matrices
#' @description
#' This function aligns the rows of a list of discrete-valued source matrices.
#' Each matrix in the list is reordered so that its rows best match the
#' corresponding rows of the first matrix in the list, which is used as the
#' reference matrix.
#'
#' @param Slist A list of discrete-valued source matrices. All matrices should
#'   have the same dimensions, with rows corresponding to latent sources and
#'   columns corresponding to variables, node pairs, or features.
#'
#' @return A list of discrete-valued source matrices with the same length as
#'   \code{Slist}. The first matrix is unchanged, and the rows of all subsequent
#'   matrices are reordered to match the rows of the first matrix.
#' @export
#'
#' @importFrom gtools permutations
#'
match_rows = function(Slist){
  S0 = Slist[[1]]
  Slist_out = list()
  Slist_out[[1]] = S0
  for(i in 2:length(Slist)){
    Slist_out[[i]] = Slist[[i]][compute.match.idx.s.tra2(Slist[[i]],S0),]
  }
  return(Slist_out)
}



#' @export
match_source = function(S_true, S_est){
  return(S_est[compute.match.idx.s.tra2(S_est,S_true),])
}

