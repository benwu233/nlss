#' @title Heatmaps for Networks
#' @description The function plots heatmaps for single or multiple network(s).
#'
#' @param S a matrix with each row representing a vectorized weighted adjacency matrix.
#' @param lim a 2-dimentional vector specifying the limits for the data,
#' @param community a vector represents the community each node belongs to.
#' @param color colorbar for the heatmap,
#' @param legend type of the legend, it can be NULL (no legend), "FC" (a symmetric, functional-connectivity-type legend) or "SC" (a asymmetric, structural-connectivity-type legend).
#' @param path the path that heatmaps will be stored at, ended with a "/" or "\".
#' @param filename names for heatmaps.
#'
#' @return
#' @export
#' @import gplots
#' @import grDevices
#' @importFrom plotrix gradient.rect
#' @importFrom graphics text
#'
#' @examples
heatmap.net = function(S,lim = c(min(S),max(S)),community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),color = bluered(100),legend = NULL, portion=FALSE, thr = 0.5,path,filename){

  sidecolor = rep("#b7b7b7",length(community))
  colsep0 = NULL
  target = "#8c8c8c"
  tmp = "#b7b7b7"
  for(i in 2:length(community)){
    if(community[i]!=community[i-1]){
      sidecolor[i] = target
      target = tmp
      tmp = sidecolor[i]
      colsep0 = c(colsep0,i-1)
    }
    else{
      sidecolor[i] = sidecolor[i-1]
    }
  }
  if(class(S)!="matrix"){
    S = t( as.matrix(S) )
  }

  if(portion==TRUE){
    if(legend=="FC"){
      minS = -1
      maxS = 1
    }
    if(legend=="SC"){
      minS = 0
      maxS = 1
    }

  }
  else{
    minS = lim[1]
    maxS = lim[2]
  }

  h = (maxS - minS)/length(color)

  for(i in 1:nrow(S)){
    ad_S = vec_mat(as.numeric(S[i,]) )
    L = nrow(ad_S)

    if(portion==TRUE){
      portion_S = matrix(0, nrow = L, ncol = L)
      group0 = unique(community)
      for(j in group0){
        ind1 = which(community==j)
        for(k in group0){
          ind2 = which(community==k)
          tmp1 = (sum( ad_S[ind1,ind2] > 0 ) + (j==k)*length(ind1) ) / length(ind1)/length(ind2)
          tmp2 = (sum( ad_S[ind1,ind2] <0 )  ) / length(ind1)/length(ind2)
          if( (tmp1 > thr) + (tmp2 > thr) > 0 ){
            if(tmp1 > tmp2){
              portion_S[ind1, ind2] = tmp1
            }
            else{
              portion_S[ind1, ind2] = -tmp2
            }
          }
        }
      }
      ad_S = portion_S
    }

    while (dev.cur()>1) dev.off()
    png(paste0(path,filename,"_",i,".png"),width = 1000, height = 850)
    heatmap.2(ad_S,Rowv = NULL,Colv = NULL,trace="none", symbreaks = TRUE,
              labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor, ColSideColors = sidecolor,
              cexRow=1,cexCol=1,margins=c(1,1),col=color, breaks=seq(minS,maxS,h),
              colsep = colsep0, rowsep = colsep0, sepcolor="black",
              key = F, lhei = c(0.75,3.5), lwid = c(1.5,3.5))
    if(is.null(legend)==FALSE){
      if(legend=="FC"){
        gradient.rect(0.065,0.2,0.08,0.75,nslices = length(color),border = T,gradient = "y",
                      col = color)
        text(x = rep(0.06,3), y = c(0.2, 0.475, 0.75),adj = 1,cex = 5,
             labels = c(minS,"0",maxS))
      }
      if(legend=="SC"){
        gradient.rect(0.065,0.2,0.08,0.75,nslices = length(color),border = T,gradient = "y",
                      col = color)
        text(x = rep(0.06,3), y = c(0.2, 0.75),adj = 1,cex = 5,
             labels = c("0",maxS))
      }
    }

    while (dev.cur()>1) dev.off()
  }
}


#' @export
heatmap.net.pp = function(S,S_sample,thr = 0.0,p = 0.95,community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),color = bluered(100), path,filename){

  sidecolor = rep("#b7b7b7",length(community))
  colsep0 = NULL
  target = "#8c8c8c"
  tmp = "#b7b7b7"
  for(i in 2:length(community)){
    if(community[i]!=community[i-1]){
      sidecolor[i] = target
      target = tmp
      tmp = sidecolor[i]
      colsep0 = c(colsep0,i-1)
    }
    else{
      sidecolor[i] = sidecolor[i-1]
    }
  }
  if(class(S)!="matrix"){
    S = t( as.matrix(S) )
  }


  minS = -1
  maxS = 1

  h = (maxS - minS)/length(color)

  n_sample = dim(S_sample)[3]

  for(i in 1:nrow(S)){
    ad_S = vec_mat(as.numeric(S[i,]) )
    L = nrow(ad_S)
    portion_S = matrix(0, nrow = L, ncol = L)
    group0 = unique(community)


    S_list = list()
    for(l in 1:n_sample){
      S_list[[l]] = vec_mat(as.numeric(S_sample[i,,l]) )
    }

    for(j in group0){
      ind1 = which(community==j)
      for(k in group0){
        ind2 = which(community==k)
        if(j <=k){
          tmp1 = (sum( ad_S[ind1,ind2] > 0 ) + (j==k)*length(ind1) ) / length(ind1)/length(ind2)
          tmp2 = (sum( ad_S[ind1,ind2] <0 )  ) / length(ind1)/length(ind2)
          if( (tmp1 > thr) + (tmp2 > thr) > 0 ){
            if(tmp1 > tmp2){
              class0 = 1
              tmp3 = tmp1
            }
            else{
              class0 = -1
              tmp3 = tmp2
            }

            tmp_ad_S = ad_S[ind1,ind2]
            sum0 = sum(tmp_ad_S==class0)
            count0 = 0
            for(l in 1:n_sample){
              tmp_S = S_list[[l]]
              tmp_S1 = tmp_S[ind1,ind2]
              if( sum(tmp_S1[tmp_ad_S==class0] == class0) >= (p* sum0) ){
                count0 = count0 + 1
              }
            }
            if(class0<0){
              count0 = -count0
            }

            portion_S[ind1,ind2] = count0 / n_sample #* tmp3
          }
        }
        else{
          portion_S[ind1,ind2] = portion_S[ind2,ind1]
        }
      }
    }
    ad_S = portion_S

    while (dev.cur()>1) dev.off()
    png(paste0(path,filename,"_",i,".png"),width = 1000, height = 850)
    heatmap.2(ad_S,Rowv = NULL,Colv = NULL,trace="none", symbreaks = TRUE,
              labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor, ColSideColors = sidecolor,
              cexRow=1,cexCol=1,margins=c(1,1),col=color, breaks=seq(minS,maxS,h),
              colsep = colsep0, rowsep = colsep0, sepcolor="black",
              key = F, lhei = c(0.75,3.5), lwid = c(1.5,3.5))

    gradient.rect(0.065,0.2,0.08,0.75,nslices = length(color),border = T,gradient = "y",
                  col = color)
    text(x = rep(0.06,3), y = c(0.2, 0.475, 0.75),adj = 1,cex = 5,
         labels = c(minS,"0",maxS))

    while (dev.cur()>1) dev.off()
  }
}



