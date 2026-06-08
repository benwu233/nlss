#' @title Heatmaps for Network Matrices
#' @description
#' This function plots heatmaps for one or multiple networks. Each row of
#' \code{S} is treated as a vectorized weighted adjacency matrix and is converted
#' back to a symmetric matrix before plotting. If multiple rows are provided,
#' the corresponding network heatmaps are displayed side by side.
#'
#' @param S A matrix or vector representing one or multiple vectorized weighted
#'   adjacency matrices. If \code{S} is a matrix, each row corresponds to one
#'   vectorized network. If \code{S} is a vector, it is treated as a single
#'   vectorized network. The length of each vectorized network should be
#'   \eqn{d(d-1)/2}, where \eqn{d} is the number of nodes.
#' @param lim A numeric vector of length two specifying the lower and upper
#'   limits of the color scale. The default is \code{c(min(S), max(S))}.
#' @param community A vector of community labels for network nodes. Its length
#'   should be equal to the number of nodes in the network. The community labels
#'   are used to add side colors and separation lines in the heatmap. The default
#'   assumes all nodes belong to the same community.
#' @param color A vector of colors used for the heatmap color scale. The default
#'   is \code{bluered(100)}.
#' @return This function is called for its side effect of producing heatmap
#'   plots. It does not return a user-level object.
#' @export
#' @import gplots
#' @import grDevices
#' @import gridGraphics
#' @import gridExtra
#' @import grid
#' @importFrom plotrix gradient.rect
#' @importFrom graphics text
#'
heatmap.net = function(S,lim = c(min(S),max(S)),
                       community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),
                       color = bluered(100)){

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
  }


  gl = lapply(1:nrow(S), function(i){
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
  grid.arrange(grobs=gl, ncol=nrow(S), clip=TRUE)
}


# heatmap.net2 = function(S,lim = c(min(S),max(S)),
#                        community = rep(1,(1 + sqrt(1+8*ncol(S))) / 2),
#                        color = bluered(100),
#                        ncol = nrow(S),
#                        filename, width, height){
#
#   sidecolor = rep("#b7b7b7",length(community))
#   colsep0 = NULL
#   target = "#8c8c8c"
#   tmp = "#b7b7b7"
#   for(i in 2:length(community)){
#     if(community[i]!=community[i-1]){
#       sidecolor[i] = target
#       target = tmp
#       tmp = sidecolor[i]
#       colsep0 = c(colsep0,i-1)
#     }
#     else{
#       sidecolor[i] = sidecolor[i-1]
#     }
#   }
#   if(class(S)[1]!="matrix"){
#     S = t( as.matrix(S) )
#   }
#
#   minS = lim[1]
#   maxS = lim[2]
#
#   h = (maxS - minS)/length(color)
#
#   ad_S = list()
#   for(i in 1:nrow(S)){
#     ad_S[[i]] = vec_mat(as.numeric(S[i,]) )
#     L = nrow(ad_S)
#     png(filename = paste(filename,"_",i,".png"), width = width, height = height)
#     heatmap.2(ad_S[[i]],Rowv = FALSE,Colv = FALSE,trace="none", symbreaks = TRUE,
#               labRow=NA,labCol=NA,dendrogram="none",RowSideColors = sidecolor,
#               ColSideColors = sidecolor,
#               cexRow=1,cexCol=1,colRow="white",colCol="white",
#               margins=c(0.2,0.2),col=color, breaks=seq(minS,maxS,h),
#               colsep = colsep0, rowsep = colsep0, sepcolor="black",
#               key = FALSE, lhei = c(0.05,3.5), lwid = c(0.05,3.5))
#     dev.off()
#   }
#
# }

