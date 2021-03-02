#' @title Simulation of NLSS
#' @description This function simulates observations using NLSS with pre-specified latent sources.
#' The number of latent sources is three, with one representing within-community connections and two
#' representing cross-comminity connections.
#'
#' @param n_node number of nodes.
#' @param n sample size.
#' @param alpha relative signal strengh.
#' @param beta relative noise strengh.
#'
#' @return a R list containing the following terms:
#' \describe{
#'   \item{X}{the observed matrix.}
#'   \item{S}{the latent source matrix.}
#'   \item{A}{the mixing coefficent matrix.}
#'   \item{community}{the community index.}
#' }
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

  sim_X = Simu_DLSS(S+1,A,K) - 1

  return(list(X = sim_X, S = S, A = A, community = modu) )

}

#' @export
simnlss = function(A,S,K){
  sim_X = Simu_DLSS(S+1,A,K) - 1
  return(sim_X)
}

#' @export
sim_ICA = function(n_node = 50, n = 40, modu, sd= 0.1){

  p = (n_node-1)*n_node/2
  q = 4
  S = matrix(0,nrow = q , ncol = p)

  A = matrix(rnorm(n*(q-1),0,1),nrow = n, ncol = q-1)
  for(i in 1:n){
    A[i,] = A[i,] / sqrt(sum(A[i,]^2) )
  }
  A = cbind(A,rep(1,n))

  #A = rdirichlet(n,c(0.5*rep(1,q-1),0.1))
  n_modu = length(table(modu))

  adM_S = matrix(0, nrow = n_node, ncol = n_node)
  for(i in 1:n_node){
    for(j in 1:n_node){
      if(modu[i] == modu[j]){
        adM_S[i,j] = 1
      }
    }
  }

  S[1,] = vec_mat(adM_S) #+ rnorm(p,mean=0,sd=sd)



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


  S[2,] = vec_mat(adM_S) #+ rnorm(p,mean=0,sd=sd)


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
  S[3,] = vec_mat(adM_S) #+ rnorm(p,mean=0,sd=sd)

  S[4,] = rnorm(p,mean=0,sd=sd) #+ sample(c(0,1),p,replace=TRUE,prob = c(0.9,0.1))

  X = A%*%S

  return(list(X=X, A=A, S=S[-nrow(S),],community = modu) )
}
