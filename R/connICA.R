#' connICA
#'
#' @param data Data matrix
#' @param q a
#'
#' @return
#'
#' @importFrom ica icafast
#' @importFrom stats quantile
#' @importFrom stats cor
#' @importFrom stats lm
#'
#' @examples
#' @export
connICA = function(data,q){
  Ss = matrix(0,nrow = ncol(data), ncol = q*100)
  for(i in 1:100){
    print(paste0("fastICA:",i,"/100." ))
    tmp = icafast(t(data),q,center =FALSE)
    Ss[,(q*(i-1)+1):(q*i)] = tmp$S
  }

  res = cor(Ss,Ss)
  res_tag = abs(res)>0.75
  for(i in 1:(q*100)){
    for(j in 1:100){
      tmp0 = res[i,(q*(j-1)+1):(q*j)]
      if( sum(res_tag[i,(q*(j-1)+1):(q*j)]) > 1) {
        res_tag[i,(q*(j-1)+1):(q*j)] = rep(0,q)
        res_tag[i,(q*(j-1)+which.max(abs(tmp0)))] = 1
      }
    }
  }

  group_list = list()
  tag = 1
  for(i in 1:(q*100)){
    if(sum(res_tag[i,])>75){
      group_list[[tag]] = which(res_tag[i,]==1)
      tag = tag + 1
    }
  }
  group_list1 = unique(group_list)
  tag = rep(1,length(group_list1))

  group_list2 = list()
  list_tag = 1
  while(sum(tag)>0){
    print(paste0("First pick-up, ", sum(tag), " remain." ))
    for(i in 1:length(tag)){
      if(tag[i]==1){
        curr_group = group_list1[[i]]
        tag[i] = 0
        break
      }
    }
    if(sum(tag)==0){
      group_list2[[list_tag]] = curr_group
      break
    }
    for(j in (i+1):length(tag)){
      if(tag[j]==1){
        ints = intersect(curr_group,group_list1[[j]])
        if(length(ints) >= 75){
          curr_group = ints
          tag[j] = 0
        }
      }
    }
    group_list2[[list_tag]] = curr_group
    list_tag = list_tag +1
  }

  tag = rep(1,length(group_list2))

  group_list3 = list()
  list_tag = 1
  while(sum(tag)>0){
    print(paste0("Second pick-up, ", sum(tag), " remain." ))
    for(i in 1:length(tag)){
      if(tag[i]==1){
        discards = 0
        curr_group = group_list2[[i]]
        tag[i] = 0
        break
      }
    }
    if(sum(tag)==0){
      group_list3[[list_tag]] = curr_group
      break
    }
    for(j in (i+1):length(tag)){
      if(tag[j]==1){
        ints = intersect(curr_group,group_list2[[j]])
        if(length(ints) >= 25){
          if(length(group_list2[[j]])<length(curr_group)) {
            tag[j] = 0
          }
          else{
            discards = 1
            break
          }
        }
      }
    }
    if(discards==0){
      group_list3[[list_tag]] = curr_group
      list_tag = list_tag +1
    }

  }

  group_list4 = list()
  list_tag = 1
  for(i in 1:length(group_list3)){
    print(paste0("checking group ",i,"/",length(group_list3),"...be patient.." ))
    Ss_0 = group_list3[[i]]
    if(sum(res_tag[Ss_0,Ss_0])==(length(Ss_0)^2) ){
      group_list4[[list_tag]] = group_list3[[i]]
      list_tag = list_tag + 1
    }
    else{
      while(sum(res_tag[Ss_0,Ss_0])<=(length(Ss_0)^2)){
        Ss_0 = Ss_0[-which.min(apply(res_tag[Ss_0,Ss_0],1,sum))]
        if(length(Ss_0)<75){
          break
        }
      }
      if(length(Ss_0)>=75) {
        group_list4[[list_tag]] = Ss_0
        list_tag = list_tag + 1
      }
    }
  }

  print(length(group_list3))
  group_mean = matrix(0,nrow = ncol(data), ncol = length(group_list4))
  for(i in 1:ncol(group_mean)){
    base0 = group_list4[[i]][1]
    groupi = Ss[,base0]
    for(j in 2:length(group_list4[[i]])){
      if(res[base0,group_list4[[i]][j]]>0) {
        groupi = cbind(groupi, Ss[,group_list4[[i]][j]])
      }
      else{
        groupi = cbind(groupi, -Ss[,group_list4[[i]][j]])
      }
    }
    group_mean[,i] = apply(groupi,1,mean)

  }

  A = matrix(0,nrow = nrow(data), ncol = ncol(group_mean))
  for(i in 1:nrow(data)){
    lm0 = lm(data[i,]~0+group_mean)
    A[i,] = lm0$coefficients
  }

  return(list(S = t(group_mean), A = A))
}
