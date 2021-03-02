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
heatmap.net.new = function(S,lim = c(min(S),max(S)),community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),color = bluered(100),legend = NULL, portion=FALSE, thr = 0.5,path,filename){

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
        gradient.rect(0.15,0.35,0.165,0.75,nslices = length(color),border = T,gradient = "y",
                      col = color)
        text(x = rep(0.14,3), y = c(0.35, 0.55, 0.75),adj = 1,cex = 5,
             labels = c(minS,"0",maxS))
      }
      if(legend=="SC"){
        gradient.rect(0.15,0.35,0.165,0.75,nslices = length(color),border = T,gradient = "y",
                      col = color)
        text(x = rep(0.14,3), y = c(0.35, 0.75),adj = 1,cex = 5,
             labels = c("0",maxS))
      }
    }

    while (dev.cur()>1) dev.off()
  }
}

#' @export
heatmap.net.new.nolg = function(S,lim = c(min(S),max(S)),community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),color = bluered(100),path,filename){

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


  minS = lim[1]
  maxS = lim[2]

  h = (maxS - minS)/length(color)

  for(i in 1:nrow(S)){
    ad_S = vec_mat(as.numeric(S[i,]) )
    L = nrow(ad_S)

    while (dev.cur()>1) dev.off()
    png(paste0(path,filename,"_",i,".png"),width = 1000, height = 850)
    heatmap.2(ad_S,Rowv = NULL,Colv = NULL,trace="none", symbreaks = TRUE,
              labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor, ColSideColors = sidecolor,
              cexRow=1,cexCol=1,margins=c(1,1),col=color, breaks=seq(minS,maxS,h),
              colsep = colsep0, rowsep = colsep0, sepcolor="black",
              key = F, lhei = c(0.75,3.5), lwid = c(1.5,3.55))

    while (dev.cur()>1) dev.off()
  }
}


#' @export
heatmap.net.new.bar = function(S,lim = c(min(S),max(S)),community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),color = bluered(100),path,filename){

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


  minS = lim[1]
  maxS = lim[2]

  h = (maxS - minS)/length(color)

  for(i in 1:nrow(S)){
    ad_S = vec_mat(as.numeric(S[i,]) )
    L = nrow(ad_S)

    while (dev.cur()>1) dev.off()
    png(paste0(path,filename,"_",i,".png"),width = 800, height = 800)
    heatmap.2(ad_S,Rowv = NULL,Colv = NULL,trace="none", symbreaks = TRUE,
              labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor, ColSideColors = sidecolor,
              cexRow=1,cexCol=1,margins=c(1,1),col=color, breaks=seq(minS,maxS,h),
              colsep = colsep0, rowsep = colsep0, sepcolor="black",
              key = F, lhei = c(1.00,3.00), lwid = c(1.00,3.00))
    gradient.rect(0.05,0.15,0.13,0.75,nslices = length(color),border = T,gradient = "y",
                  col = color)

    while (dev.cur()>1) dev.off()
  }
}
