#------------------------------------------------------------------------------#
#-------------------------- Common Functions ----------------------------------#
#------------------------------------------------------------------------------#

# INPUT: 
# n <- no. of samples
# dim <- dimension of the data
# OUTPUT:
# n many samples from the uniform distribution

# Function for generating samples under null

X.gen <- function(n,dim){
  if (dim == 1){
    X.samp <- as.matrix(replicate(n, runif(dim,0,1)))
  }
  else{
    X.samp <- t(replicate(n, runif(dim,0,1)))
  }
  return (X.samp)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Function for generating samples from perturbed uniform distribution
# INPUT:
# n <- sample size
# dum <- dimension of samples
# num_perturb <- number of perturbations

Y.gen <- function(n, dim, num_perturb, sob_smooth = 1, c_d = 2.7, 
                  theta_seed = 100){
  theta_seed <- sample(1:(1000*n), 1)
  samp_seed <- sample(1:(1000*n), 1)
  Y.samp <- f_theta_sampler(theta_seed, samp_seed, n, num_perturb, sob_smooth,
                            c_d, dim)
  
  return (Y.samp)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Function for computing MMD^2 between two data vectors
# INPUTS:
# k <- given kernel
# X,Y <- given dataset
# OUTPUTS:
# MMD.out <- value of MMD_{u}^{2} for given kernel k

compute.MMD <- function(X, Y, k){
  
  # Compute K(Xi,Xj), K(Yi,Yj) and K(Xi,Yj)
  k.X <- kernlab::kernelMatrix(k, X)
  k.Y <- kernlab::kernelMatrix(k, Y)
  k.XY <- kernlab::kernelMatrix(k, X, Y)
  
  # Compute MMD_{u}^{2} statistic based on above values
  MMD.out <- mean(k.X[row(k.X)!=col(k.X)] + k.Y[row(k.Y)!=col(k.Y)] -
                    2*k.XY[row(k.XY)!=col(k.XY)])
  
  return (MMD.out)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for computing vector of MMD_{u}^{2} values for a given finite number
# of kernels
# INPUT:
# kernel.vec <- a list of given kernels
# X,Y <- given dataset
# OUTPUT:
# MMD.vec <- vector of MMD_{u}^{2} for list of given kernels

compute.MMD.vec <- function(X, Y, kernel.vec){
  MMD.vec <- rep(NA, length(kernel.vec))
  # looping over each kernel
  for (i in 1:length(kernel.vec)){
    # MMD_{u}^2 value for the i^{th} kernel
    MMD.val <- compute.MMD(X, Y, kernel.vec[[i]])
    MMD.vec[i] <- MMD.val
  }
  return (MMD.vec)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for estimating the asymptotic covariance matrix under H0 using only X
# INPUTS
# n <- sample size
# x <- data
# k.vec <- list of kernels
# OUTPUT
# The estimated covariance matrix

est.cov <- function(n,x,k.vec){
  # required for using trace function
  require(psych)
  # length of the list of kernels
  k.len <- length(k.vec)
  # Centering matrix required for centering the kernel matrix
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  # list of kernel matrices indexed by the list of kernels
  kvec.mat <- vector("list", k.len)
  for (i in 1:k.len){
    # Computing the centered kernel matrix
    kvec.mat[[i]] <- C%*%kernlab::kernelMatrix(k.vec[[i]], x)%*%C
  }
  
  # estimated covariance matrix
  cov.mat.est <- matrix(0, nrow = k.len, ncol = k.len)
  for (i in 1:k.len){
    for (j in 1:k.len){
      cov.mat.est[i,j] <- (8/(n^2))*tr(kvec.mat[[i]]%*%kvec.mat[[j]])
    }
  }
  
  # adding small error for ensuring invertibility
  return (cov.mat.est + (10^-5)*min(diag(cov.mat.est))*diag(1, k.len,k.len))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for choosing median bandwidth given a data
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
  # parameter is taken as exp(-sigma*d(x,y)), but
  # according to median heurestic paper the "sigma" 
  # is standard, hence I have this relation.
  
  return (sigma.hat)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------- Single Kernel Functions ----------------------------#
#------------------------------------------------------------------------------#

# Function for returning the cut-off for test using single kernel. Here we use a
# formulation slightly different from gretton, instead of using the eigenvalues,
# we equivalently use the quadratic form.
# INPUTS: 
# n <- no. of sample of X
# x <- the data X
# k <- the kernel to use
# OUTPUT: 
# The estimated cutoff under H0.

single.H0.cutoff <- function(n, x, k, n.iter = 1000){
  # required for using trace function
  require(psych)
  # Centering matrix
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  # Centered and scaled kernel matrix
  k.mat <- (1/n)*(C%*%kernlab::kernelMatrix(k, x)%*%C)
  
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for simulating test using single kernel
# INPUTS
# n <- number of original samples;
# d <- dimension of the Gaussian data
# p <- number of perturbations
# sob_smooth <- sobolev smoothness parameter
# c_d <- multiplier for perturbed density
# theta_seed <- seed for generating theta
# kernel.choice <- choice of kernel used
# Output
# power of test using single kernel

Single.MMD <- function(n, d, p, sob_smooth = 1, c_d = 2.7, 
                       theta_seed = 1, kernel.choice = "GAUSS", n.iter = 1000){
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, sob_smooth, c_d, 
                                              theta_seed)})
  count <- foreach(k=1:n.iter, .combine= 'c',
                   .export = ls(envir=globalenv())) %dopar% {
    
                     
    cat (paste(k,"\n"), file = "log.txt", append = TRUE)
    X <- X.list[[k]]
    Y <- Y.list[[k]]
    
    # choosing rbf/laplace kernel bandwidth from data using median heurestic
    sigma.med <- med.bandwidth(X, Y)
    
    # Fixing kernel according to above choice
    if (kernel.choice == "GAUSS"){
      kernel <- kernlab::rbfdot(sigma = sigma.med)
    }
    else if (kernel.choice == "LAP"){
      kernel <- kernlab::laplacedot(sigma = sqrt(sigma.med))
    }
    
    # Threshold for 0.05 level test
    MMD.threshold <- single.H0.cutoff(n, X, kernel, n.iter)
    
    # n*MMD_{u}^{2} value for the given sample
    MMD.val <- n*compute.MMD(X,Y,kernel)
    
    (MMD.val > MMD.threshold)
    
    
  }
  print ("Single Kernel: Done")
  return (sum(count)/n.iter)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------- Multiple Kernel Functions --------------------------#
#------------------------------------------------------------------------------#

# function for choosing range of bandwidth using min-max heurestic
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
  
  return (2/quantiles.norm)    
  # From kernlab package documention, the "sigma" parameter is taken as 
  # exp(-sigma*d(x,y)), but according to median heurestic paper the "sigma" 
  # is standard, hence this formulation.
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for choosing range based on 2^l rule
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
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for computing function of MMD_{u}^{2} vector for a list of kernels
# INPUTS
# x <- a vector of values
# param <- some other parameter to be used
# OUTPUT
# value of the function applied on the vector x

multi.func <- function(x, param){
  # Now, param is set to be the cov matrix of asymptotic null distribution. 
  # Hence we consider the following Mahalanobis type statistic. Observe that 
  # under null the asymptotic mean is 0.
  
  out <- t(x)%*%param%*%x
  #out <- max(x)
  
  return (out)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for finding H0 cutoff for multiple Kernel test.
# INPUT:
# n <- number of samples from distribution of X.
# x <- the oberseved data coming from p
# k.vec <- list of kernels to use
# OUTPUT:
# The upper alpha cutoff under null.


multi.k.approx.stat <- function(k.mat, u.mat){
  
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  return (test.stat)
}

multi.H0.cutoff <- function(n, x, k.vec, n.iter = 1000){
  # required for using trace function
  require(psych)
  # length of the list of kernels
  k.len <- length(k.vec)
  
  # Centering matrix
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  # list of kernel matrices indexed by the list of kernels
  kvec.mat <- vector("list", k.len)
  for (i in 1:k.len){
    # Centered kernel matrices scaled by 1/n
    kvec.mat[[i]] <- (1/n)*(C%*%kernlab::kernelMatrix(k.vec[[i]], x)%*%C)
  }
  # Estimated covariance matrix
  invcov.mat.est <- spdinv(est.cov(n, x, k.vec))
  
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  test.kernel.mat <- sapply(kvec.mat, multi.k.approx.stat, u.mat = u.mat)
  test.stat <- apply(test.kernel.mat, 1, multi.func, param = invcov.mat.est)
  
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for providing list of kernels according to user choice
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

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Function for simulating test using multiple kernel
# INPUTS
# n <- number of original samples;
# d <- dimension of the data
# p <- number of perturbations
# sob_smooth <- sobolev smoothness parameter
# c_d <- multiplier for perturbed density
# theta_seed <- seed for generating theta
# kernel.choice <- choice of kernel used

# OUTPUT
# power of test using single kernel
Multi.MMD <- function(n, d, p, sob_smooth = 1, c_d = 2.7, 
                      theta_seed = 1, kernel.choice = "GEXP",
                      n.iter = 1000){
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, sob_smooth, c_d, 
                                              theta_seed)})
  
  count <- foreach(i=1:n.iter,.inorder = T,.combine = 'c',
                   .export = ls(envir=globalenv()), .packages = "Rfast") %dopar% {
    
    cat (paste(i,"\n"), file = "log.txt", append = TRUE)
                     
    X <- X.list[[i]]
    Y <- Y.list[[i]]
    
    # Kernels (sigmas) to use 
    kernel.vec <- k.choice(X, Y, kernel.choice)
    
    
    # Threshold for 0.05 level test using above functional (estimated using above
    # samples)
    MMD.func.threshold <- multi.H0.cutoff(n, X, kernel.vec, n.iter)
    
    # Value of functional of n*MMD_{u}^{2} vector for above list of kernels
    MMD.samp.val <- n*compute.MMD.vec(X, Y, kernel.vec)
    inv.cov.samp <- Rfast::spdinv(est.cov(n,X,kernel.vec))
    MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
    
    
    MMD.samp.func > MMD.func.threshold
    
  }
  print ("Multiple Kernel: Done")
  return (sum(count)/n.iter)
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#---------------------------- Power plot function -----------------------------#
#------------------------------------------------------------------------------#
# Function for comparing power over dimensions
# INPUTS:
# n <- no of data points from both distributions
# d <- dimension of the data
# p.seq <- vector of number of perturbations
# sob_smooth <- sobolev smoothness parameter
# c_d <- multiplier for perturbed density
# theta_seed <- seed for generating theta
# kernel.choice <- choice of kernels under single (first two coordinates) 
# and multi kernel setup (last three coordinate)
# OUTPUT:
# A data frame having powers of single and multiple test under various
# perturbations
power.d <- function(n, d, p.seq,sob_smooth,c_d,
                    kernel.choice,theta_Seed = 1,  n.iter = 500){
  
  writeLines(c(""), "log.txt")
  
  # redefining perturbations vector for ease
  p <- p.seq
  out.compare <- c()
  for (k in 1:length(p)){
    
    #--------------------------------------------------------------------------#
    library(Rfast)
    #--------------------------------------------------------------------------#
    # Estimating power under single kernel test

    out.row.col1 <- Single.MMD(n, d, p[k],sob_smooth = 1, c_d = 2.7, 
                               theta_seed = 1, kernel.choice[1], n.iter)
    cat(paste("Single Kernel-1 in iteration",out.row.col1," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    out.row.col2 <- Single.MMD(n, d, p[k], sob_smooth = 1, c_d = 2.7, 
                               theta_seed = 1, kernel.choice[2], n.iter)
    cat(paste("Single Kernel-2 in iteration",out.row.col2," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    # Estimating power under multiple kernel test
    
    out.row.col3 <- Multi.MMD(n, d, p[k],sob_smooth = 1, c_d = 2.7, 
                              theta_seed = 1,kernel.choice[3], n.iter)
    cat(paste("Multiple Kernel-1 in iteration",out.row.col3," ",k,"\n"),
        file="log.txt", append=TRUE)
    start <- Sys.time()
    out.row.col4 <- Multi.MMD(n, d, p[k],sob_smooth = 1, c_d = 2.7, 
                              theta_seed = 1,kernel.choice[4], n.iter)
    Sys.time()-start
    cat(paste("Multiple Kernel-2 in iteration",out.row.col4," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    out.row.col5 <- Multi.MMD(n, d, p[k],sob_smooth = 1, c_d = 2.7, 
                              theta_seed = 1,kernel.choice[5], n.iter)
    cat(paste("Multiple Kernel-3 in iteration",out.row.col5," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    # Concatenating all the outputs
    out.row <- c(p[k],out.row.col1, out.row.col2, out.row.col3, out.row.col4,
                 out.row.col5)
    
    out.compare <- rbind(out.compare, out.row)
    
  }
  
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Number of Perturb","Single Kernel-1","Single Kernel-2",
                             "Multiple Kernel-1", "Multiple Kernel-2",
                             "Multiple Kernel-3")
  
  return (out.compare)
}
