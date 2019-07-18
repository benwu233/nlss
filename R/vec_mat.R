#' @title transformation between a symmetric matrix and a vector
#' @description \code{vec_mat} is a small function that vectorize a symmetric
#' matrix (only keep the upper triangle) or transform a vector to a symmetric
#' matrix
#'
#' @param x a symmetric matrix or a vector
#'
#' @return a vector or a symmetric matrix
#' @export
#'
#' @examples
#' x = matrix(rnorm(36),nrow = 6)
#' vec_mat(x)
vec_mat = function(x){

  if(class(x)=="matrix"){
    n = nrow(x)*(nrow(x) -1)/2
    out = rep(NA,n)
    tag = 1
    for(i in 1:(nrow(x)-1) ){
      for(j in (i+1):ncol(x)){
        out[tag] = x[i,j]
        tag = tag + 1
      }
    }
  }
  else if(class(x)=="numeric"){
    n = (1 + sqrt(1+8*length(x))) / 2
    out = matrix(NA,ncol = n, nrow = n)
    tag = 1
    for(i in 1:(n-1) ){
      for(j in (i+1):n){
        out[i,j] = x[tag]
        tag = tag + 1
      }
    }

    diag(out) = 0
    for(i in 2:n){
      for(j in 1:(i-1) ){
        out[i,j] = out[j,i]
      }
    }
  }

  return(out)
}
