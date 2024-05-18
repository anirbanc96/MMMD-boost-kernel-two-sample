################################################################################
#################### Functions for computing MMD values ########################

########################### Compute MMD_{u}^{2} ################################

compute.MMD <- function(X, Y, k){
  
  ##############################################################################
  # input: X,Y <- given data matrix
  #        k <- considerd kernel
  # output: MMD.out <- the value of MMD_{u}^{2} between X and Y
  ##############################################################################
  
  k.X <- kernlab::kernelMatrix(k, X)
  k.Y <- kernlab::kernelMatrix(k, Y)
  k.XY <- kernlab::kernelMatrix(k, X, Y)
  
  MMD.out <- mean(k.X[row(k.X)!=col(k.X)] + k.Y[row(k.Y)!=col(k.Y)] -
                    2*k.XY[row(k.XY)!=col(k.XY)])
  
  return (MMD.out)
}

################################################################################

############## Compute MMD_{u}^{2} vector for multiple kernels #################

compute.MMD.vec <- function(X, Y, kernel.vec){
  
  ##############################################################################
  # input: X,Y <- given data matrix
  #        k <- considerd kernel list
  # output: MMD.vec <- the value of MMD_{u}^{2} between X and Y for list k
  ##############################################################################
  
  MMD.vec <- rep(NA, length(kernel.vec))
  
  for (i in 1:length(kernel.vec)){
    
    MMD.val <- compute.MMD(X, Y, kernel.vec[[i]])
    MMD.vec[i] <- MMD.val
  }
  return (MMD.vec)
}

################################################################################

##################### Estimate the covariance matrix ###########################

est.cov <- function(n,x,k.vec){
  
  ##############################################################################
  # input: n <- number of samples
  #        x <- data under null distribution
  #        k.vec <- list of kernels considered
  # output: estimated covariance matrix with small additive error
  ##############################################################################
  
  require(psych)
  k.len <- length(k.vec)
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  kvec.mat <- vector("list", k.len)
  for (i in 1:k.len){
    kvec.mat[[i]] <- C%*%kernlab::kernelMatrix(k.vec[[i]], x)%*%C
  }
  
  cov.mat.est <- matrix(0, nrow = k.len, ncol = k.len)
  for (i in 1:k.len){
    for (j in 1:k.len){
      cov.mat.est[i,j] <- (8/(n^2))*tr(kvec.mat[[i]]%*%kvec.mat[[j]])
    }
  }
  
  return (cov.mat.est + (10^-5)*min(diag(cov.mat.est))*diag(1, k.len,k.len))
}

################################################################################

########################### Median Bandwidth ###################################

med.bandwidth <- function(X, Y){
  
  ##############################################################################
  # input: X,Y <- given datasets
  # output: sigma.hat <- squared median bandwidth
  ##############################################################################
  
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  Z <- rbind(X,Y)
  
  nu.med <- median(dist(Z)^2)
  
  sigma.hat <- 1/nu.med
  ###################################################
  # From kernlab package documention, the "sigma"   #
  # parameter is taken as exp(-sigma*d(x,y)), but   #
  # according to median heurestic paper the "sigma" #
  # is standard, hence I have this relation.        #
  ###################################################
  return (sigma.hat)
}

################################################################################
################# Choice of bandwidth using min-max heurestic ##################

min.max.band <- function(X, Y){
  
  ##############################################################################
  # input: X,Y <- input datasets
  # output: chosen bandwidths
  ##############################################################################
  
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  Z <- rbind(X,Y)
  norm.vals <- dist(Z)^2
  quantiles.norm <- quantile(norm.vals, probs = c(0.95,0.05))
  
  return (1/quantiles.norm)    
}

################################################################################

##################### Choice of bandwidth using 2^l rule #######################

expo.band <- function(X,Y, l0 = -3, l1 = 3){
  
  ##############################################################################
  # input: X,Y <- input datasets
  #        l0 <- lower bound of range
  #        l1 <- upper bound of range
  # output: a vector of bandwidths
  ##############################################################################
  med.band <- med.bandwidth(X,Y)
  
  band <- c();i <- l0
  while(i<=l1){
    band <- c(band,(2^i)*med.band)
    i <- i + 1
  }
  return(band)
}

################################################################################
################## Choice of kernels for multi kernel test #####################

k.choice <- function(X,Y, kernel.choice){
  
  ##############################################################################
  # input: X,Y <- given datasets
  #        kernel.choice <- string for choice of kernels
  # output: list of kernels according to given choice
  ##############################################################################
  
  # List of kernels using min-max bandwidth
  if (kernel.choice == "MINMAX"){
    sigma.len <- 10
    sigma.bd <- min.max.band(X,Y)
    sigma.vec <- seq(sigma.bd[1],sigma.bd[2], length.out = sigma.len)
    
    kernel.vec <- c()
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec, kernlab::rbfdot(sigma = sigma.vec[i]))
    }
  }
  
  # List of kernels using exponential band method
  else if (kernel.choice == "GEXP"){
    
    l0 <- -2; l1 <- 2
    sigma.vec <- expo.band(X,Y, l0, l1)
    sigma.len <- length(sigma.vec)
    
    kernel.vec <- c()
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec, kernlab::rbfdot(sigma = sigma.vec[i]))
    }
  }
  # Mixture of RBF and LAPLACE kernels , where bandwidth is selected using 
  # exponential method
  else if (kernel.choice == "MIXED"){
    
    l0 <- -1; l1 <- 1
    sigma.vec <- expo.band(X,Y, l0, l1)
    sigma.len <- length(sigma.vec)
    
    kernel.vec <- c()
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec, kernlab::rbfdot(sigma = sigma.vec[i]))
    }
    
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec,
                      kernlab::laplacedot(sigma = sqrt(sigma.vec[i])))
    }
  }
  # Laplace kernels using exponential band method
  else if (kernel.choice == "LAP"){
    
    l0 <- -2; l1 <- 2
    sigma.vec <- expo.band(X,Y, l0, l1)
    sigma.len <- length(sigma.vec)
    
    kernel.vec <- c()
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec,
                      kernlab::laplacedot(sigma = sqrt(sigma.vec[i])))
    }
  }
  return (kernel.vec)
}

################################################################################
