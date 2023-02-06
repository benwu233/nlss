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
#' @import gridGraphics
#' @import gridExtra
#' @import grid
#' @importFrom plotrix gradient.rect
#' @importFrom graphics text
#'
#' @examples
heatmap.net = function(S,lim = c(min(S),max(S)),
                       community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),
                       color = bluered(100),
                       ncol = nrow(S)){

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
  if(class(S)[1]!="matrix"){
    S = t( as.matrix(S) )
  }

  minS = lim[1]
  maxS = lim[2]

  h = (maxS - minS)/length(color)

  ad_S = list()
  for(i in 1:nrow(S)){
    ad_S[[i]] = vec_mat(as.numeric(S[i,]) )
    L = nrow(ad_S)
  }


  gl = lapply(1:3, function(i){
    heatmap.2(ad_S[[i]],Rowv = FALSE,Colv = FALSE,trace="none", symbreaks = TRUE,
              labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor,
              ColSideColors = sidecolor,
              cexRow=1,cexCol=1,colRow="white",colCol="white",
              margins=c(0.2,0.2),col=color, breaks=seq(minS,maxS,h),
              colsep = colsep0, rowsep = colsep0, sepcolor="black",
              key = FALSE, lhei = c(0.05,3.5), lwid = c(0.05,3.5))
    grid.echo()
    grid.grab()
  })
  grid.newpage()
  grid.arrange(grobs=gl, ncol=ncol, clip=TRUE)
}

