#' @title Row-Wise Thresholding by Normal Tail Probabilities
#'
#' @description
#' This function thresholds each row of a numeric matrix using two-sided normal
#' tail probabilities. For each row, the standard deviation of that row is used
#' as the scale parameter, and entries with p-values smaller than \code{a} are
#' marked as \code{TRUE}.
#'
#' @param X A numeric matrix whose rows are thresholded separately.
#' @param a A numeric value specifying the significance threshold for the
#'   two-sided normal tail probability. The default is \code{0.05}.
#'
#' @return A logical matrix with the same dimensions as \code{X}. Entries are
#'   \code{TRUE} if the corresponding value in \code{X} has a two-sided normal
#'   tail probability smaller than \code{a}, using the row-wise standard
#'   deviation as the scale parameter, and \code{FALSE} otherwise.
#'
#' @export
thr_byrow= function(X,a=0.05){
  out = matrix(0,nrow = nrow(X),ncol = ncol(X))
  sd_vec = rep(0, ncol(X))
  for(i in 1:nrow(X)){
    p = pnorm(abs(X[i,]),sd=sd(X[i,]),lower.tail = FALSE)*2
    out[i,] = (t(p)<a)
  }
  return(out)
}
