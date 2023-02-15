################################################################################
########################### Common Functions ###################################
################################################################################
# mean under H0 <- mu0,mean under H1 <- mu1; 
# Covariance Matrix under H0 <- Sigma0, 
# Covariance Matrix under H1 <- Sigma1 
# are used as global variables and hence are not passed as function arguments


################### Function for generating samples under null #################
# INPUT: 
# n <- no. of samples
# dim <- dimension of the data
# p <- probability of mixture
# OUTPUT:
# n many samples from the distribution under H0

X.gen <- function(n, dim, p){
  
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
Y.gen <- function(n, dim, p){
  
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
########## Function for computing MMD_u^2 between two data vectors #############
# INPUTS:
# k <- given kernel (kernlab object)
# X,Y <- given dataset
# OUTPUTS:
# MMD.out <- value of MMD_{u}^{2} for given kernel k

compute.MMD <- function(X, Y, k){
  
  # Compute K(Xi,Xj), K(Yi,Yj) and K(Xi,Yj) matrices
  k.X <- kernlab::kernelMatrix(k, X)
  k.Y <- kernlab::kernelMatrix(k, Y)
  k.XY <- kernlab::kernelMatrix(k, X, Y)
  
  # Compute MMD_{u}^{2} statistic based on above values
  MMD.out <- mean(k.X[row(k.X)!=col(k.X)] + k.Y[row(k.Y)!=col(k.Y)] -
                    2*k.XY[row(k.XY)!=col(k.XY)])
  
  return (MMD.out)
}

################################################################################
# Function for computing vector of MMD_{u}^{2} values for a given finite number
# of kernels
################################################################################
# INPUT:
# kernel.vec <- a list of given kernels (vector of kernlab objects)
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
################################################################################
# Function for estimating the asymptotic covariance matrix under H0 using only 
# data coming from null distribution
################################################################################
# INPUTS
# n <- sample size
# x <- data
# k.vec <- list of kernels (list of kernlab objects)
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

################################################################################
########################## Single Kernel Functions #############################
################################################################################

####### Function for returning the cut-off for test using single kernel ########

# INPUTS: 
# n <- no. of sample of X, i.e. samples from p
# x <- the data X
# k <- the kernel to use (kernlab object)
# OUTPUT: 
# The estimated cutoff under H0.

single.H0.cutoff <- function(n, x, k, n.iter){
  # required for using trace function
  require(psych)
  # Centering matrix
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  # Centered and scaled kernel matrix
  k.mat <- (1/n)*(C%*%kernlab::kernelMatrix(k, x)%*%C)
  
  # generating all the gaussian random vectors for repeated generation
  # of samples from the approximating distribution
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  # computing n.iter many samples of approximating distribution
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  # computing upper 95% quantile of approximating distribution
  H0.thresh <- quantile(test.stat, probs = 0.95)
  
  return (H0.thresh)
}

# Observe that here we do not use the eigenvalue approach from Gretton (2014)
# instead we use the equivalent quadratic form. This approach is exactly what
# outlined in Gretton (2014) as using property of Gaussian distribution
# the quadratic form using eigenvalues can equivalently written as a 
# quadratic form using the centered kernel matrix

################################################################################
############## Function for simulating test using single kernel ################
# INPUTS
# n <- number of original samples;
# d <- dimension of the data
# p <- the probability of mixing.
# kernel.choice <- choice of kernel considered ("GAUSS" or "LAP")
# n.iter <- number of iterations to be done.

# Output
# proportion of rejections

Single.MMD <- function(n, d, p, kernel.choice, n.iter){
  
  # variable for couting total number of rejections
  count <- 0
  
  for (i in 1:n.iter){
    
    # generating samples from H0 and H1
    X <- X.gen(n, d, p)
    Y <- Y.gen(n, d, p)
    
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
    
    # Threshold for 0.05 level test
    MMD.threshold <- single.H0.cutoff(n, X, kernel, n.iter)
    
    # n*MMD_{u}^{2} value for the given sample
    MMD.val <- n*compute.MMD(X,Y,kernel)
    
    # Adding rejection to previous count of rejections
    count <- count + (MMD.val > MMD.threshold)
    #print (c(i,count))
  }
  return (count/n.iter)
}
################################################################################
################################################################################

################################################################################
########################## Multiple Kernel Functions ###########################
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

expo.band <- function(X,Y, l0, l1){
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
# param <- a matrix of appropriate dimensions
# OUTPUT
# value of the function applied on the vector x

multi.func <- function(x, param){
  
  out <- t(x)%*%param%*%x
  
  return (out)
}
################################################################################
############ Function for finding H0 cutoff for multiple Kernel test ###########
# INPUT:
# n <- number of samples from distribution of X.
# x <- the oberseved data coming from p
# k.vec <- list of kernels to use (list of kernlab objects)
# invcov <- the estimated inverse covariance matrix
# n.iter <- number of iterations for finding threshold
# OUTPUT:
# The upper alpha cutoff under null.

####### Helper function for computing value of approximating statistic given a #
# kernel and the independently generated gaussian variables ####################
multi.k.approx.stat <- function(k.mat, u.mat){
  
  # computing n.iter many samples of approximating distribution for a single
  # kernel
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  return (test.stat)
}

#---------------------######################################-------------------#

multi.H0.cutoff <- function(n, x, k.vec, invcov, n.iter){
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
  # Estimated covariance matrix (which is passed as an argument)
  invcov.mat.est <- invcov
  
  # generating all the gaussian random vectors for repeated generation
  # of samples from the approximating distribution
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  # getting n.iter many resamples from the approximating distribution
  # first we generate a matrix having n.iter many rows, where the i^th
  # row contains the random vector generated from the approximating distribution
  # using the i^th gaussian sample
  test.kernel.mat <- sapply(kvec.mat, multi.k.approx.stat, u.mat = u.mat)
  
  # from the above generated matrix, we apply the mahalanobis distance type
  # function to each row with the estimated covariance matrix
  test.stat <- apply(test.kernel.mat, 1, multi.func, param = invcov.mat.est)
  
  # computing upper 95% quantile of approximating distribution
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}

################################################################################
###### Function for providing list of kernels according to user choice #########
# INPUT: 
# X,Y <- observed data (used for finding bandwidth)
# kernel.choice <- choice of kernel ("MINMAX", "GEXP", "MIXED", "LAP")
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
############ Function for simulating test using multiple kernel ################
# INPUTS
# n <- number of original samples
# d <- dimension of the data
# p <- mixing probability
# kernel.choice <- a string representing the chosen kernel method 
#                  ("MINMAX", "GEXP", "MIXED", "LAP")
# n.iter <- number of iterations for finding threshold

# OUTPUT
# proportion of iterations rejected
Multi.MMD <- function(n, d, p, kernel.choice, n.iter){
  
  # variable for couting total number of rejections
  count <- 0
  
  # n.iter sized lists containing the generated data to be used in each
  # iteration
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p)})
  
  for (i in 1:n.iter){
    
    # Calling i^th samples from above pre-generated data list
    X <- X.list[[i]]
    Y <- Y.list[[i]]
    
    # Kernels to use 
    kernel.vec <- k.choice(X, Y, kernel.choice)
    
    
    # Threshold for 0.05 level test using above function
    
    # estimated, corrected for singularity inverse covariance matrix
    inv.cov.samp <- spdinv(est.cov(n,X,kernel.vec))
    
    # Storing the threshold coming from above function
    MMD.func.threshold <- multi.H0.cutoff(n, X, kernel.vec, 
                                          inv.cov.samp, n.iter)
    
    # Value of n*MMD_{u}^{2} vector for above list of kernels
    MMD.samp.val <- n*compute.MMD.vec(X, Y, kernel.vec)
    # Value of the mahalanobish type functional applied to above sample vector
    MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
    
    # Adding rejection to previous count of rejections
    count <- count + (MMD.samp.func > MMD.func.threshold)
  }
  return (count/n.iter)
}


################################################################################
################################################################################

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
# n.iter <- number of iterations for estimating power.
# OUTPUT:
# A data frame having powers of single and multiple test under various
# dimensions
power.d <- function(n, sigma.param, sigma.mult, mu.param, d.seq, p,
                    kernel.choice, n.iter = 500){
  
  # redefining dimension vector for ease
  d <- d.seq
  
  out.compare <- foreach(k=1:length(d), .combine=rbind, .export = ls(envir=globalenv())) %dopar% {
    
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
    # cov matrix under H0
    Sigma0 <- matrix(0, nrow = d[k], ncol = d[k])
    for(i in 1:d[k]){
      for(j in 1:d[k]){
        Sigma0[i,j] = sigma.param^(abs(i-j))
      }
    }
    
    ##### Uncomment below line and comment above section for diagonal Cov ######
    
    #Sigma0 <- diag(sigma.param, d[k], d[k])
    
    # cov matrix under H1
    Sigma1 <- sigma.mult*Sigma0
    
    # mean vector under H0
    mu0 <- rep(0, d[k])
    # mean vector under H1
    mu1 <- rep(mu.param, d[k])
    
    # Estimating power under single kernel test
    
    out.row.col1 <- Single.MMD(n, d[k], p,kernel.choice[1], n.iter)
    
    out.row.col2 <- Single.MMD(n, d[k], p,kernel.choice[2], n.iter)
    
    # Estimating power under multiple kernel test
    
    out.row.col3 <- Multi.MMD(n, d[k], p,kernel.choice[3], n.iter)
    
    out.row.col4 <- Multi.MMD(n, d[k], p,kernel.choice[4], n.iter)
    
    out.row.col5 <- Multi.MMD(n, d[k], p,kernel.choice[5], n.iter)
    
    # Concatenating all the outputs
    out.row <- c(d[k],out.row.col1, out.row.col2, out.row.col3, out.row.col4,
                 out.row.col5)
  }
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Dimension","Single Kernel-1","Single Kernel-2",
                             "Multiple Kernel-1", "Multiple Kernel-2",
                             "Multiple Kernel-3")
  
  return (out.compare)
}
