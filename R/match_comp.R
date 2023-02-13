#' @title Match discrete value matrices by their rows
#' @description The function reorders the rows of matrices to match them to rows of the
#' first matrix in the list.
#'
#' @param Slist list of discrete value source matrices with same dimensions.
#'
#' @return A list of discrete value source matrices with rows matched to the rows of the first matrix.
#' @export
#'
#' @examples
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

