#' @title Vectorization of a List of Symmetric Matrices
#' @description
#' This function vectorizes a list of symmetric matrices by extracting the
#' lower triangular entries of each matrix and stacking the resulting vectors
#' row by row. Optionally, each vectorized matrix can be thresholded by its own
#' empirical quantile and further binarized.
#'
#' @param mat_list A list of symmetric matrices. All matrices should have the
#'   same dimensions.
#' @param q A numeric value between 0 and 1 specifying the quantile threshold
#'   applied to each vectorized matrix. Entries smaller than or equal to the
#'   \code{q}-th empirical quantile are set to zero. The default is \code{0},
#'   in which case no thresholding is applied.
#' @param binarize A logical value indicating whether the output should be
#'   binarized. If \code{TRUE}, nonzero entries are returned as \code{TRUE} and
#'   zero entries as \code{FALSE}. The default is \code{FALSE}.
#'
#' @return A matrix in which each row corresponds to one matrix in
#'   \code{mat_list}, and each column corresponds to one lower triangular entry
#'   of the original symmetric matrices. If \code{binarize = TRUE}, a logical
#'   matrix is returned.
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
