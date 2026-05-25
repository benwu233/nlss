#' @title Generate latent sources
#' @description This function simulates generate three latent source networks with eight communities.
#' @param n_node number of nodes.
#' @return a R list containing the following terms:
#' \describe{
#'   \item{S}{the latent source matrix with each row representing a source network (vectorized adjacency matrix).}
#'   \item{community}{the community index for each node.}
#' }
#' @export
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
gen_sources = function(n_node = 50){

  p = (n_node-1)*n_node/2
  q = 3
  S = matrix(0,nrow = q, ncol = p)

  modu = rep(0,n_node)
  for( i in 1:n_node){
    if(i <= (n_node/10*1.5) ){
      modu[i] = 1
    }
    else if(i <= (n_node/10*3.5)){
      modu[i] = 2
    }
    else if(i <= (n_node/10*5)){
      modu[i] = 3
    }
    else if(i <= (n_node/10*5.5)){
      modu[i] = 4
    }
    else if(i <= (n_node/10*6.5)){
      modu[i] = 5
    }
    else if(i <= (n_node/10*8.5)){
      modu[i] = 6
    }
    else if(i <= (n_node/10*9)){
      modu[i] = 7
    }
    else{
      modu[i] = 8
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
      if( ((modu[i] == 3)  +(modu[j] == 3) ) > 0){
        adM_S[i,j] = 1
      }
    }
  }

  S[2,] = vec_mat(adM_S)

  adM_S = matrix(0, nrow = n_node, ncol = n_node)
  for(i in 1:n_node){
    for(j in 1:n_node){
      if( ((modu[j] == 7) + (modu[i] == 7) )> 0){
        adM_S[i,j] = 1
      }
    }
  }
  S[3,] = vec_mat(adM_S)

  return(list(S=S, community = modu) )
}

#' @title Simulate data from a NLSS model adding Gaussian noise
#' @description This function simulates observations using NLSS with given latent sources.
#' @param n sample size
#' @param alpha_0 concentration parameter (a three dimensional vector) for latent sources.
#' @param alpha_1 concentration parameter (a scale variable) for noise.
#' @param sd0 standard deviation for the Gaussian noise
#' @param S the latent source matrix
#' @return a R list containing the following terms:
#' \describe{
#'   \item{Xc}{the data matrix (with Gaussian noise)}
#'   \item{X}{the data matrix (without Gaussian noise)}
#'   \item{A}{the mixing coefficient matrix}
#' }
#' @export
sim_NLSS = function(n = 50, alpha_0 = c(0.5,0.5,0.5), alpha_1 = 0.1,
                    sd0 = 0.3, S){

  q = nrow(S)
  p = ncol(S)
  A = rdirichlet(n,c(alpha_0,alpha_1))

  sim_X = simNLSS(S,A,unique(as.numeric(S)))

  noise0 = matrix(rnorm(n*p), nrow = n, ncol =p)
  sim_X1 = sim_X + sd0*noise0

  return(list(Xc = sim_X1, X = sim_X, A = A) )
}

#' @title Simulate data from a linear mixture with Gaussian noise
#' @description This function simulates observations using a linear mixture with given latent sources.
#' @param n sample size
#' @param sd standard deviation for the Gaussian mixing coefficient
#' @param sd0 standard deviation for the Gaussian noise
#' @param S the latent source matrix
#' @return a R list containing the following terms:
#' \describe{
#'   \item{X}{the data matrix}
#'   \item{A}{the mixing coefficient matrix}
#'   \item{noise}{the noise matrix}
#' }
#' @export
sim_ICA = function(n = 50, sd = c(0.5,0.5,0.5),
                  sd0= 0.1, S){

  p = ncol(S)

  q = nrow(S)
  A = matrix(0,nrow = n, ncol = q)

  for(j in 1:q){
    A[,j] = rnorm(n,0,sd[j])
  }

  X = A%*%S

  noise0 =  matrix(rnorm(n*p,mean = 0, sd = sd0), nrow = n, ncol = p)

  X = X + noise0

  return(list(X=X, A=A, noise=noise0) )
}




