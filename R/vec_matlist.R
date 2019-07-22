#' @title Vectrization of a list of symmetric matrices
#' @description Vectorize a list of symmetric matrices and stack them together.
#' The outcome can be thresholded and binarized with optional parameters.
#'
#' @param mat_list a list of matrices.
#' @param q a quantile serving as a threshold, any value smaller than the quantile will be set to zero.
#' @param binarize logical. Binarize the outcome or not.
#'
#' @return a matrix with each row representing a matrix from the input mat_list.
#' @export
#'
#' @examples
vec_matlist = function(mat_list, q = 0.0, binarize = FALSE){

  out = NULL

  if(q>0){
    for(i in 1:length(mat_list)){

      tmp = vec_mat(mat_list[[i]])
      threshold = quantile(tmp,q)

      out = rbind( out,t( tmp*(tmp>threshold) )  )
    }
  }
  else{
    for(i in 1:length(mat_list)){

      tmp = vec_mat(mat_list[[i]])

      out = rbind( out,tmp )
    }
  }

  rownames(out) <- names(mat_list)

  if(binarize){
    return(out!=0)
  }
  else{
    return(out)
  }
}
