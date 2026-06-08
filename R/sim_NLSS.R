#' @title Generate Latent Source Networks
#' @description
#' This function generates three latent source networks with a predefined
#' eight-community structure. The generated source networks include a
#' within-community source and two hub-like sources associated with selected
#' communities.
#' @param n_node An integer specifying the number of nodes in the network.
#'   The default is \code{50}.
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{S}}{A latent source matrix, where each row represents one
#'   source network in vectorized adjacency matrix form.}
#'   \item{\code{community}}{A vector of community labels for the network nodes.}
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

#' @title Simulate Data from an NLSS Model with Gaussian Noise
#' @description
#' This function simulates observations from the Network Latent Source
#' Separation (NLSS) model with given latent source networks. Mixing coefficients
#' are generated from a Dirichlet distribution, and optional Gaussian noise is
#' added to the noiseless NLSS observations.
#' @param n An integer specifying the sample size. The default is \code{50}.
#' @param alpha_0 A numeric vector of Dirichlet concentration parameters for the
#'   latent source components. Its length should be equal to the number of
#'   latent sources, namely \code{nrow(S)}. The default is
#'   \code{c(0.5, 0.5, 0.5)}.
#' @param alpha_1 A numeric value specifying the Dirichlet concentration
#'   parameter for the background or noise component. The default is \code{0.1}.
#' @param sd0 A numeric value specifying the standard deviation of the additive
#'   Gaussian noise. The default is \code{0.3}.
#' @param S A latent source matrix, where each row represents one source network
#'   in vectorized adjacency matrix form.
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{Xc}}{The simulated data matrix with additive Gaussian noise.}
#'   \item{\code{X}}{The noiseless simulated data matrix generated from the
#'   NLSS model.}
#'   \item{\code{A}}{The simulated mixing coefficient matrix.}
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

#' @title Simulate Data from a Linear Mixture with Gaussian Noise
#'
#' @description
#' This function simulates observations from a linear mixture model with given
#' latent source networks. The mixing coefficients are generated independently
#' from zero-mean Gaussian distributions, and Gaussian noise is added to the
#' resulting linear mixture.
#'
#' @param n An integer specifying the sample size. The default is \code{50}.
#' @param sd A numeric vector specifying the standard deviations of the Gaussian
#'   mixing coefficients for different latent sources. Its length should be
#'   equal to the number of latent sources, namely \code{nrow(S)}. The default is
#'   \code{c(0.5, 0.5, 0.5)}.
#' @param sd0 A numeric value specifying the standard deviation of the additive
#'   Gaussian noise. The default is \code{0.1}.
#' @param S A latent source matrix, where each row represents one source network
#'   in vectorized adjacency matrix form.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{X}}{The simulated data matrix from the noisy linear mixture
#'   model.}
#'   \item{\code{A}}{The simulated mixing coefficient matrix.}
#'   \item{\code{noise}}{The additive Gaussian noise matrix.}
#' }
#'
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




