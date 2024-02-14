################### Function for generating samples under null #################
# INPUT: 
# n <- no. of samples
# dim <- dimension of the data
# p <- probability of mixture
# OUTPUT:
# n many samples from the distribution under H0

X.gen <- function(n, dim, p, mu0, mu1, 
                  Sigma0, Sigma1){
  
  # generating sample from gaussian denoted by mixture probability = 0
  if (p == 0){
    X.samp <- MASS::mvrnorm(n, mu = mu0, Sigma = Sigma0)
  }
  # generating sample from t (df=10) or multivariate Laplace distribution
  # denoted by mixture probability = 1. 
  else if (p == 1){
    X.samp <- LaplacesDemon::rmvt(n, mu = mu0, S = Sigma0, df = 10)
    
    ############################################################################
    # Uncomment line below and comment line above for generating from Laplace
    
    #X.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
    ############################################################################
  }
  else{
    
    # vector of choices for mixture distributions
    choice.vec <- rbinom(n, 1, p)
    
    # generating the samples from Gaussian distribution
    gauss.samp <- MASS::mvrnorm(n, mu = mu0, Sigma = Sigma0)
    
    # generating the samples from t (df = 10) or Laplace for mixture
    other.samp <- LaplacesDemon::rmvt(n, mu = mu0, S = Sigma0, df = 10)
    ############################################################################
    # Uncomment line below and comment line above for generating from Laplace
    
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
    ############################################################################
    X.samp <- (1-choice.vec) * gauss.samp + choice.vec * other.samp
  }
  return (X.samp)
}

################################################################################
############### Function for generating samples under alternative ##############
# INPUT: 
# n <- no. of samples
# dim <- dimension of the data
# p <- probability of mixture
# OUTPUT:
# n many samples from the distribution under H1
Y.gen <- function(n, dim, p, mu0, mu1, 
                  Sigma0, Sigma1){
  
  # generating sample from gaussian denoted by mixture probability = 0
  if (p == 0){
    Y.samp <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma1)
  }
  # generating sample from t (df=10) or multivariate Laplace distribution
  # denoted by mixture probability = 1. 
  else if (p == 1){
    Y.samp <- LaplacesDemon::rmvt(n, mu = mu1, S = Sigma1, df = 10)
    ############################################################################
    # Uncomment line below and comment line above for generating from Laplace
    
    #Y.samp  <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
  }
  else{
    
    # vector of choices for mixture distributions
    choice.vec <- rbinom(n, 1, p)
    
    # generating the samples from Gaussian distribution
    gauss.samp <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma1)
    
    # generating the samples from t (df = 10) or Laplace for mixture
    other.samp <- LaplacesDemon::rmvt(n, mu = mu1, S = Sigma1, df = 10)
    ############################################################################
    # Uncomment line below and comment line above for generating from Laplace
    
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
    
    Y.samp <- (1-choice.vec) * gauss.samp + choice.vec * other.samp
  }
  return (Y.samp)
}

################################################################################
############### The Function $L_{K}$ for a given kernel K ######################
# INPUT:
# d <- dimension of the data
# z1 <- first datavector in the form of [x,y] in a rowbind format
# z2 <- second datavector in the form of [x,y] in a rowbind format
# k <- chosen kernel
# OUTPUT:
# the value of the function $L_{k}$.
L_k <- function(d, z1, z2, k){
  
  return (k(z1[1:d], z2[1:d]) + k(z1[(d+1):(2*d)], z2[(d+1):(2*d)])
          - k(z1[1:d],z2[(d+1):(2*d)]) - k(z2[1:d], z1[(d+1):(2*d)]))
  
}

################################################################################
#################### Function to compute $MMD_{\ell}$ ##########################
# INPUT:
# n <- number of samples
# d <- dimension of the data
# Z1 <- first data matrix having rows in the form [x,y] 
# Z2 <- second data matrix having rows in the form [x,y]
# kernel <- the chosen kernel
# OUTPUT:
# the value of $MMD_{\ell}$ based on above data
compute.MMD_l <- function(n,d,Z1, Z2, kernel){
  
  L_k.vec <- apply(cbind(Z1,Z2), 1,
                   function(x){L_k(d, x[1:(2*d)],x[(2*d+1):(4*d)], kernel)})
  
  return (mean(L_k.vec))
  
}

################################################################################
###################

compute.MMD_l.vec <- function(n,d,Z1, Z2, kernel.vec){
  
  r <- length(kernel.vec)
  MMD_l.vec <- rep(NA, r)
  for (i in 1:r){
    
    MMD_l.vec[i] <- compute.MMD_l(n,d,Z1, Z2, kernel.vec[[i]])
    
  }
  
  return (MMD_l.vec)
}

est.cov.ab <- function(n, d, Z1, Z2, ka, kb){
  
  L1 <- apply(cbind(Z1, Z2), 1,
              function(x){L_k(d, x[1:(2*d)],x[(2*d+1):(4*d)], ka)})
  L2 <- apply(cbind(Z1, Z2), 1,
              function(x){L_k(d, x[1:(2*d)],x[(2*d+1):(4*d)], kb)})
  
  sigma.ab <- cov(L1, L2)
  
  return (2*sigma.ab)
  
}

est.cov <- function(n,d, Z1, Z2, k.vec){
  
  r <- length(k.vec)
  cov.mat <- matrix(NA, r, r)
  for (i in 1:r){
    for (j in 1:r){
      
      cov.mat[i,j] <- est.cov.ab(n,d, Z1, Z2, k.vec[[i]], k.vec[[j]])
      
    }
  }
  
  return (cov.mat + (10^-5)*min(diag(cov.mat))*diag(1, r,r))
}


################################################################################
############ Function for choosing median bandwidth given data #################
# INPUTS:
# X,Y <- given dataset
# OUTPUTS:
# sigma.hat <- bandwidth (for rbf/laplace in kernlab) using median bandwidth

med.bandwidth <- function(X, Y){
  
  # making X and Y matrices for generality
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  Z <- rbind(X,Y)
  
  # median of the row distances
  nu.med <- median(dist(Z)^2)
  
  # rbfdot kernel bandwidth
  sigma.hat <- 1/nu.med        
  # From kernlab package documention, the "sigma"
  # parameter for rbf is taken as exp(-sigma*d(x,y)^2), but
  # according to median heurestic the "sigma" 
  # is standard, hence I have this relation.
  
  return (sigma.hat)
}

# Note: The above function actually returns the square of median bandwidth. This
# is the primary reason why we use square root of output of above function for
# Laplace kernels. 

################################################################################
################################################################################

###### Function for choosing range of bandwidth using min-max heurestic ########
# INPUTS
# X,Y <- given dataset
# OUTPUT
# upper and lower bound of range of bandwidth using min-max heurestic

min.max.band <- function(X, Y){
  
  # making X and Y matrices for generality
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  Z <- rbind(X,Y)
  
  # distance between rows of above data matrix
  norm.vals <- dist(Z)^2
  quantiles.norm <- quantile(norm.vals, probs = c(0.95,0.05))
  
  return (1/quantiles.norm)    
  # From kernlab package documention, the "sigma" parameter for rbf is taken as 
  # exp(-sigma*d(x,y)^2), but according to median heurestic paper the "sigma" 
  # is standard, hence this formulation.
}

################################################################################
############## Function for choosing range based on 2^l rule ###################
# INPUTS
# X,Y <- given dataset
# l0, l1 <- upper and lower bound of integer range
# OUTPUT
# a vector of bandwidths

expo.band <- function(X,Y, l0 = -3, l1 = 3){
  # Computing median bandwidth
  med.band <- med.bandwidth(X,Y)
  
  band <- c();i <- l0
  # loop for choosing bandwidths as i ranges over all integers between l0 and l1
  while(i<=l1){
    band <- c(band,(2^i)*med.band)
    i <- i + 1
  }
  return(band)
}

################################################################################
# Function for computing function of MMD_{u}^{2} vector for a list of kernels ##
################################################################################
# INPUTS
# x <- a vector of values
# param <- the estimated and singularity corrected covariance matrix
# OUTPUT
# value of the function applied on the vector x

multi.func <- function(x, param){
  
  out <- t(x)%*%param%*%x
  
  return (out)
}
################################################################################
################################################################################
###### Function for providing list of kernels according to user choice #########
# INPUT: 
# X,Y <- observed data (used for finding bandwidth)
# kernel.choice <- choice of kernel (between rbf laplace or mixed)
# OUTPUT:
# A list containing kernels using a pre-specified bandwidth selection method
k.choice <- function(X,Y, kernel.choice){
  
  # List of kernels using min-max bandwidth
  if (kernel.choice == "MINMAX"){
    sigma.len <- 10
    sigma.bd <- min.max.band(X,Y)
    sigma.vec <- seq(sigma.bd[1],sigma.bd[2], length.out = sigma.len)
    
    # Defining and Storing list of kernels
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
    
    # Defining and Storing list of kernels
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
    
    # Defining and Storing list of kernels
    kernel.vec <- c()
    for(i in 1:sigma.len){
      kernel.vec <- c(kernel.vec,
                      kernlab::laplacedot(sigma = sqrt(sigma.vec[i])))
    }
  }
  
  # Return the list of chosen kernels.
  return (kernel.vec)
}

################################################################################
########### Function for simulating test based on single kernel ################


Single.MMD <- function(n, d, p, kernel.choice, n.iter, mu0, mu1, 
                       Sigma0, Sigma1){
  
  # variable for couting total number of rejections
  count <- 0
  cat(" ", file = "log.txt")
  # n.iter sized lists containing the generated data to be used in each
  # iteration
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p,mu0, mu1, 
                                              Sigma0, Sigma1)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, mu0, mu1, 
                                              Sigma0, Sigma1)})
  
  count.vec <- foreach(i=1:n.iter, .combine=cbind, .export = ls(envir=globalenv())) %dopar% {
    
    cat(c("Started Iteration", i, "\n"), file = "log.txt", append = T)
    library(Rfast)
    
    X <- X.list[[i]]
    Y <- Y.list[[i]]
    
    Z <- cbind(X, Y)
    Z1 <- Z[seq(1, n, 2),]
    Z2 <- Z[seq(2, n, 2),]
    
    # Kernels to use 
    #kernel.vec <- k.choice(X, Y, kernel.choice)
    
    # choosing rbf/laplace kernel bandwidth from data using median heurestic
    sigma.med <- med.bandwidth(X, Y)
    
    # Fixing kernel according to above choice
    if (kernel.choice == "GAUSS"){
      kernel <- kernlab::rbfdot(sigma = sigma.med)
    }
    else if (kernel.choice == "LAP"){
      # Note that we take square-root of the output of median bandwidth function
      # This is because we output the median of squared distances, and for 
      # Laplace kernel we need the median of the distances and not the 
      # squared distances
      kernel <- kernlab::laplacedot(sigma = sqrt(sigma.med))
    }
    
    # Threshold for 0.05 level test using above function
    
    # estimated, corrected for singularity inverse covariance matrix
    #inv.cov.samp <- spdinv(est.cov(n, d, Z1,Z2,kernel.vec))
    
    inv.cov.samp <- est.cov.ab(n, d, Z1, Z2, kernel, kernel)
    
    # Storing the threshold coming from above function
    #MMD.func.threshold <- qchisq(0.95, df = length(kernel.vec))
    
    MMD.func.threshold <- qchisq(0.95, df = 1)
    
    # Value of n*MMD_{u}^{2} vector for above list of kernels
    MMD.samp.val <- sqrt(n)*compute.MMD_l(n,d,Z1, Z2, kernel)
    # Value of the mahalanobish type functional applied to above sample vector
    MMD.samp.func <- (MMD.samp.val/sqrt(inv.cov.samp))^2
    
    cat(c("Completed Single Iteration", i, "\n"), file = "log.txt", append = T)
    
    # Adding rejection to previous count of rejections
    MMD.samp.func > MMD.func.threshold
    
  }
  return (mean(count.vec))
}


################################################################################
############ Function for simulating test using multiple kernel ################
# INPUTS
# n <- number of original samples;
# d <- dimension of the data
# p <- mixing probability
# X,Y <- given dataset
# kernel.choice <- a string representing the chosen kernel method

# OUTPUT
# proportion of iterations rejected
Multi.MMD <- function(n, d, p, kernel.choice, n.iter, mu0, mu1, 
                      Sigma0, Sigma1){
  
  # variable for couting total number of rejections
  count <- 0
  cat(" ", file = "log.txt")
  # n.iter sized lists containing the generated data to be used in each
  # iteration
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p,mu0, mu1, 
                                              Sigma0, Sigma1)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, mu0, mu1, 
                                              Sigma0, Sigma1)})
  
  count.vec <- foreach(i=1:n.iter, .combine=cbind, .export = ls(envir=globalenv())) %dopar% {
    
    cat(c("Started Iteration", i, "\n"), file = "log.txt", append = T)
    library(Rfast)
    
    X <- X.list[[i]]
    Y <- Y.list[[i]]
    
    Z <- cbind(X, Y)
    Z1 <- Z[seq(1, n, 2),]
    Z2 <- Z[seq(2, n, 2),]
    
    # Kernels to use 
    kernel.vec <- k.choice(X, Y, kernel.choice)
    
    
    # Threshold for 0.05 level test using above function
    
    # estimated, corrected for singularity inverse covariance matrix
    inv.cov.samp <- spdinv(est.cov(n, d, Z1,Z2,kernel.vec))
    
    # Storing the threshold coming from above function
    MMD.func.threshold <- qchisq(0.95, df = length(kernel.vec))
    
    # Value of n*MMD_{u}^{2} vector for above list of kernels
    MMD.samp.val <- sqrt(n)*compute.MMD_l.vec(n,d,Z1, Z2, kernel.vec)
    # Value of the mahalanobish type functional applied to above sample vector
    MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
    
    cat(c("Completed Iteration", i, "\n"), file = "log.txt", append = T)
    
    # Adding rejection to previous count of rejections
    MMD.samp.func > MMD.func.threshold
    
  }
  return (mean(count.vec))
}

################################################################################
########################## Power comparison function ###########################
################################################################################
# INPUTS:
# n <- no of data points from both distributions
# sigma.param <- the parameter used for generating variance covariance matrix
# sigma.mult <- the parameter for multiplying cov matrix under H0 to get cov
# under H1
# mu.param <- the parameter for generating mean under H1
# d.seq <- the vector of dimensions
# p <- the mixing probability
# kernel.choice <- choice of kernels under single (first two coordinates) 
# and multi kernel setup (last three coordinates)
# OUTPUT:
# A data frame having powers of single and multiple test under various
# dimensions
power.d <- function(n, sigma.param, sigma.mult, mu.param, d, p,
                    kernel.choice, n.iter = 500){
  # writing a log file to keep track of progress
  out.compare <- matrix(NA, nrow = length(sigma.mult), 6)
  for (k in 1:length(sigma.mult)){
    
    #--------------------------------------------------------------------------#
    library(LaplacesDemon)
    library(Rfast)
    #--------------------------------------------------------------------------#
    
    #--------------------------------------------------------------------------#
    # Creating covariance matrices under H0 and H1.
    #################################### Note ##################################
    # We keep a provision for using the correlated as well as the diagonal
    # covariance matrices. Once should be careful to note that for the diagonal
    # case the parameter sigma.param becomes the diagonal elements of the
    # covariance matrix and hence it should be made equal to the variance and
    # not the standard deviation in that case
    ############################################################################
    # Sigma0 matrix <- cov matrix under H0
    #Sigma0 <- matrix(0, nrow = d, ncol = d)
    #for(i in 1:d){
    #  for(j in 1:d){
    #    Sigma0[i,j] = sigma.param^(abs(i-j))
    #  }
    #}
    
    ##### Uncomment below line and comment above section for diagonal Cov ######
    
    Sigma0 <- diag(sigma.param, d, d)
    
    #---------------------------######################-------------------------#
    # cov matrix under H1
    Sigma1 <- sigma.mult[k]*Sigma0
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    # mean vector under H0
    mu0 <- rep(0, d)
    # mean vector under H1
    mu1 <- rep(mu.param, d)
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    
    
    out.row.col1 <- Multi.MMD(n, d, p,kernel.choice[1], n.iter, mu0, mu1, 
                              Sigma0, Sigma1)
    
    out.row.col2 <- Multi.MMD(n, d, p,kernel.choice[2], n.iter, mu0, mu1, 
                              Sigma0, Sigma1)
    
    out.row.col3 <- Multi.MMD(n, d, p,kernel.choice[3], n.iter, mu0, mu1, 
                              Sigma0, Sigma1)
    
    out.row.col4 <- Single.MMD(n, d, p,kernel.choice[4], n.iter, mu0, mu1, 
                              Sigma0, Sigma1)
    
    out.row.col5 <- Single.MMD(n, d, p,kernel.choice[5], n.iter, mu0, mu1, 
                              Sigma0, Sigma1)
    # Concatenating all the outputs
    out.row <- c(sigma.mult[k],out.row.col1, out.row.col2,
                 out.row.col3, out.row.col4, out.row.col5)
    
    out.compare[k,] <- out.row
  }
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Signal Strength",
                             "Multiple Kernel-1", "Multiple Kernel-2",
                             "Multiple Kernel-3", "Single Kernel-1",
                             "Single Kernel-2")
  
  return (out.compare)
}

