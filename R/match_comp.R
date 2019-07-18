#' Title
#'
#' @param Slist list of q*p discrete value source matrics
#'
#' @return list of q*p discrete value source matrics with rows matched to  the first matrix
#' @export
#'
#' @examples
#' @importFrom gtools permutations
match_comp = function(Slist){
  S0 = Slist[[1]]
  Slist_out = list()
  Slist_out[[1]] = S0
  for(i in 2:length(Slist)){
    Slist_out[[i]] = Slist[[i]][compute.match.idx.s.tra2(Slist[[i]],S0),]
  }

  return(Slist_out)
}
