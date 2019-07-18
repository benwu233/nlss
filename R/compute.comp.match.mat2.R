compute.comp.match.mat2 = function(S_true,S_est){
  idx = expand.grid(1:nrow(S_true),1:nrow(S_est))


  return(matrix(sapply_pb(1:nrow(idx),
                          function(i) mean(S_true[idx[i,1],]==S_est[idx[i,2],])),
                nrow=nrow(S_true),ncol=nrow(S_est)))

}
