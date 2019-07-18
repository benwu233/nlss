#' Title
#'
#' @param n_node number of nodes
#' @param n sample size
#' @param alpha alpha
#' @param beta beta
#'
#' @return
#' @export
#'
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @examples
sim_NLSS = function(n_node = 50, n = 40, alpha = 0.5, beta = 0.1 ){

  p = (n_node-1)*n_node/2
  q = 4
  S = matrix(0,nrow = q-1, ncol = p)
  A = matrix(0,nrow = n, ncol = q)
  K = 2

  A[,1] = runif(n, min = 0.3, max = 0.5)
  A[,2] = runif(n, min = 0.15, max = 0.2)
  A[,3] = runif(n, min = 0.15, max = 0.2)
  A[,4] = 1 - A[,1] - A[,2] - A[,3]

  A = rdirichlet(n,c(alpha*rep(1,q-1),beta))

  modu = rep(0,n_node)
  for( i in 1:n_node){
    if(i <= (n_node/5*2) ){
      modu[i] = 1
    }
    else if(i <= (n_node/5*3)){
      modu[i] = 2
    }
    else{
      modu[i] = 3
    }
  }

  adM_S = matrix(0, nrow = n_node, ncol = n_node)
  for(i in 1:n_node){
    for(j in 1:n_node){
      if(modu[i] == modu[j]){
        adM_S[i,j] = 1
      }
    }
  }
  S[1,] = vec_mat(adM_S)

  adM_S = matrix(0, nrow = n_node, ncol = n_node)
  for(i in 1:n_node){
    for(j in 1:n_node){
      if( ((modu[i] == 2)+(modu[j] == 2)) > 0){
        if(modu[i]!=modu[j]){
          adM_S[i,j] = 1
        }
      }
    }
  }

  S[2,] = vec_mat(adM_S)

  adM_S = matrix(0, nrow = n_node, ncol = n_node)
  for(i in 1:n_node){
    for(j in 1:n_node){
      if( ((modu[i] == 1)*(modu[j] == 3) + (modu[i] == 3)*(modu[j] == 1)) > 0){
        if(modu[i]!=modu[j]){
          adM_S[i,j] = 1
        }
      }
    }
  }

  S[3,] = vec_mat(adM_S)

  sim_X = Simu_DLSS(S+1,A,K)

  return(list(X = sim_X, S = S, A = A, community = modu) )

}
