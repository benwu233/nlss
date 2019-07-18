compute.match.idx.s.tra2 = function(S_est,S_true){

  dist_mat = compute.comp.match.mat2(S_true,S_est)
  per = permutations(ncol(dist_mat),nrow(dist_mat))

  sum.match = function(x){
    out = 0
    for(i in 1:length(x)){
      out = out + dist_mat[i,x[i]]
    }
    return(out)
  }

  return( per[which.max(apply(per,1,sum.match)),] )
}
