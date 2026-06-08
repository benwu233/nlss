#' @title Transform Between a Symmetric Matrix and a Vector
#' @description
#' \code{vec_mat} converts a symmetric matrix to a vector by extracting its
#' lower triangular entries, or converts a vector back to a symmetric matrix by
#' filling the lower triangular part and then symmetrizing the matrix.
#'
#' @param x A symmetric matrix or a numeric vector. If \code{x} is a matrix, its
#'   lower triangular entries, excluding the diagonal, are returned as a vector.
#'   If \code{x} is a numeric vector, its length should be of the form
#'   \eqn{d(d - 1) / 2} for some integer \eqn{d}, and it is converted to a
#'   \eqn{d \times d} symmetric matrix with zero diagonal.
#'
#' @return A numeric vector or a symmetric matrix. If \code{x} is a matrix, the
#'   function returns a vector containing the lower triangular entries of
#'   \code{x}. If \code{x} is a vector, the function returns the corresponding
#'   symmetric matrix with zero diagonal.
#'
#' @examples
#' x <- matrix(rnorm(36), nrow = 6)
#' x <- (x + t(x)) / 2
#' v <- vec_mat(x)
#' X <- vec_mat(v)
#' @export
vec_mat = function(x){

  if(class(x)[1]=="matrix"){
    out = x[lower.tri(x,diag=FALSE)]
  }
  else if(class(x)[1]=="numeric"){
    n = (1 + sqrt(1+8*length(x))) / 2
    out = matrix(0,ncol = n, nrow = n)

    out[lower.tri(out,diag=FALSE)] = x

    out = out + t(out)
  }

  return(out)
}


#' @export
vec_mat0 = function(x){

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
    #out = x[upper.tri(x,diag=FALSE)]
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

    #out[upper.tri(out,diag=FALSE)] = x

   # out = out + t(out)
  }

  return(out)
}

