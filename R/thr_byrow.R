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
