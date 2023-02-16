################################################################################
########################### Common Functions ###################################
################################################################################

################################################################################
################ Function for MMD_{u}^{2}(\bm{X},\bm{Y}, k) ####################

# INPUTS:
# X,Y <- given datasets
# k <- given kernel (kernlab objects)

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
################################################################################

################################################################################
########## Function for MMD_{u}^{2}(\bm{X},\bm{Y}; k_{1},..,k_{r}) #############

# INPUT:
# X,Y <- given dataset
# kernel.vec <- a list of given kernels (kernlab objects)

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

################################################################################
################## Estimated Covariance Matrix under H0 ########################

# INPUTS
# n <- sample size
# x <- data
# k.vec <- list of kernels (kernlab objects)

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

################################################################################
############################ Median Bandwidth ##################################

# INPUTS:
# X,Y <- given dataset

# OUTPUTS:
# sigma.hat <- squared median bandwidth

med.bandwidth <- function(X, Y){
  
  # making X and Y matrices for generality
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  Z <- rbind(X,Y)
  
  # median of the row distances
  nu.med <- median(dist(Z)^2)
  
  # rbfdot kernel bandwidth
  sigma.hat <- 1/nu.med        
  
  # Note: The above formulation follows because of the usage of bandwidth in
  # kernlab package. Look at documentation of kernlab package for more details
  
  return (sigma.hat)
}

################################################################################
################################################################################


################################################################################
########################## Single Kernel Functions #############################
################################################################################

# INPUTS: 
# m <- no. of sample from H0 (here the code is done assuming balanced condition)
# x <- the data X
# k <- the kernel to use (kernlab objects)
# n.iter <- number of iterations done to estimate H0 cut-off

# OUTPUT: 
# The estimated cutoff under H0.

single.H0.cutoff <- function(m, x, k, n.iter = 1000){
  
  # required for using trace function
  require(psych)
  
  # Centering matrix
  C <- diag(1, nrow = m, ncol = m) - (1/m)*matrix(1, nrow = m, ncol = m)
  
  # Centered and scaled kernel matrix
  k.mat <- (1/m)*(C%*%kernlab::kernelMatrix(k, x)%*%C)
  
  # generating all the gaussian random vectors for repeated generation
  # of samples from the approximating distribution
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,m), 
                         Sigma = diag(2, nrow = m, ncol = m))
  
  # computing n.iter many samples of approximating distribution
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  # computing upper 95% quantile of approximating distribution
  H0.thresh <- sort(test.stat)[as.integer(0.95*n.iter)]
  
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
# m, n <- number of samples in the H0 and H1 datasets;
# resamp.size <- number of resamples considered
# kernel.choice <- choice of the kernel considered ("GAUSS", "LAP")
# n.iter <- number of iterations to be done for estimating power

# Output
# proportion of iterations rejected

Single.MMD <- function(m,n,resamp.size, kernel.choice, n.iter = 1000){
  
  count <- 0
  
  for (i in 1:n.iter){
    
    # resampled datasets from full data under H0 and H1
    
    X <- as.matrix(data.X[resamp.x[i,],])
    Y <- as.matrix(data.Y[resamp.y[i,],])
    
    # renaming the number of resamples considered
    n.samp <- resamp.size
    
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
    MMD.threshold <- single.H0.cutoff(n.samp, X, kernel, n.iter)
    
    # n*MMD_{u}^{2} value for the given sample
    MMD.val <- n.samp*compute.MMD(X,Y,kernel)
    
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
  # From kernlab package documention, the "sigma" parameter is taken as 
  # exp(-sigma*d(x,y)), but according to median heurestic paper the "sigma" 
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
# invcov <- the estimated inverse covariance matrix
# k.vec <- list of kernels to use (kernlab objects)
# OUTPUT:
# The upper alpha cutoff under null.

####### Helper function for computing value of approximating statistic given a #
# kernel and the independently generated gaussian variables ####################

multi.k.approx.stat <- function(k.mat, u.mat){
  
  # computing n.iter many samples of approximating distribution for a kernel
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  return (test.stat)
}

multi.H0.cutoff <- function(n, x, k.vec, invcov, n.iter = 1000){
  
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
  H0.thresh <- sort(test.stat)[as.integer(0.95*n.iter)]
  return (H0.thresh)
}

################################################################################
###### Function for providing list of kernels according to user choice #########
# INPUT: 
# X,Y <- observed data (used for finding bandwidth)
# kernel.choice <- choice of kernel ("MINMAX", "GEXP", "MIXED" or "LAP")
# OUTPUT:
# A list containing kernels using a pre-specified bandwidth selection method

k.choice <- function(X,Y, kernel.choice){
  
  # List of kernels using min-max bandwidth
  if (kernel.choice == "MINMAX"){
    sigma.len <- 5
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
      kernel.vec <- c(kernel.vec, kernlab::laplacedot(sigma = sqrt(sigma.vec[i])))
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
      kernel.vec <- c(kernel.vec, kernlab::laplacedot(sigma = sqrt(sigma.vec[i])))
    }
  }
  
  # Return the list of chosen kernels.
  return (kernel.vec)
}

################################################################################
############ Function for simulating test using multiple kernel ################

# INPUTS
# m, n <- number of samples in the H0 and H1 datasets;
# resamp.size <- number of resamples considered
# kernel.choice <- choice of the kernel considered 
#                  ("MINMAX", "GEXP", "MIXED" or "LAP")
# n.iter <- number of iterations to be done for estimating power

# Output
# proportion of iterations rejected

Multi.MMD <- function(m,n, resamp.size, kernel.choice,
                      n.iter = 1000){
  
  count <- 0
  
  for (i in 1:n.iter){
    
    # resampled datasets from full data under H0 and H1
    X <- as.matrix(data.X[resamp.x[i,],])
    Y <- as.matrix(data.Y[resamp.y[i,],])
    
    # Kernels (sigmas) to use 
    kernel.vec <- k.choice(X, Y, kernel.choice)
    
    
    # Threshold for 0.05 level test using above function
    
    # estimated, corrected for singularity inverse covariance matrix
    inv.cov.samp <- spdinv(est.cov(resamp.size,X,kernel.vec))
    
    # Storing the threshold coming from above function
    MMD.func.threshold <- multi.H0.cutoff(resamp.size, X, kernel.vec,
                                          inv.cov.samp, n.iter)
    
    # Value of n*MMD_{u}^{2} vector for above list of kernels
    MMD.samp.val <- resamp.size*compute.MMD.vec(X, Y, kernel.vec)
    
    # Value of the mahalanobish type functional applied to above sample vector
    MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
    
    # Adding rejection to previous count of rejections
    count <- count + (MMD.samp.func > MMD.func.threshold)
    #print (c(i,count))
  }
  return (count/n.iter)
}

################################################################################
################################################################################

################################################################################
########################## Power comparison function ###########################
################################################################################

# Function for comparing power over dimensions
# INPUTS:
# resamp <- number of resamples to be done
# set.choice <- vector of sets chosen
# kernel.choice <- choice of kernels under single (first two coordinates) 
# and multi kernel setup (last three coordinates)
# n.iter <- number of iterations to be done for estimating H0 threshold and 
# estimate power.

# OUTPUT:
# A data frame having powers of single and multiple test under various
# noise variance values
power.d <- function(resamp, set.choice,
                    kernel.choice, n.iter = 500){
  #----------------------------------------------------------------------------#
  # Libraries for parallelising
  library(foreach)
  library(doParallel)
  #----------------------------------------------------------------------------#
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  
  registerDoParallel(cl)
  #----------------------------------------------------------------------------#
  writeLines(c(""), "log.txt")
  #----------------------------------------------------------------------------#
  out.compare <- foreach(k=1:length(set.choice), .combine=rbind,
                         .export = ls(envir=globalenv())) %dopar% {
    
    #--------------------------------------------------------------------------#
    library(LaplacesDemon)
    library(Rfast)
    library(SpatialPack)
    #--------------------------------------------------------------------------#
    # choice of set for H0 and H1                       
    set.x <- set.choice[[k]][1,]; set.y <- set.choice[[k]][2,]
    #--------------------------------------------------------------------------#
    # getting indices of above set choices
    index.y <- rep(F,length(train.label)); index.x <- rep(F,length(train.label))
    for(i in 1:length(train.label)){
      if (train.label[i] %in% set.y){
        index.y[i] <- T
      }
      if (train.label[i] %in% set.x){
        index.x[i] <- T
      }
    }
    #--------------------------------------------------------------------------#
    # data under H0
    data.X <- train.x[index.x,]
    # data under H1
    data.Y <- train.x[index.y,]
    #--------------------------------------------------------------------------#
    # size of data under H0 and H1
    m.x <- dim(data.X)[1]
    n.y <- dim(data.Y)[1]
    #--------------------------------------------------------------------------#
    # indices under every resamples which is done "resamp" many times
    resamp.x <- matrix(NA, nrow = n.iter, ncol = resamp)
    resamp.y <- matrix(NA, nrow = n.iter, ncol = resamp)
    
    for(i in 1:n.iter){resamp.x[i,] <- sample(1:m.x, resamp, rep = T)}
    for(i in 1:n.iter){resamp.y[i,] <- sample(1:n.y, resamp, rep = T)}
    
    #--------------------------------------------------------------------------#
    # Creating a log file to keep track of progress
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",k,"\n"))
    #--------------------------------------------------------------------------#
    
    #--------------------------------------------------------------------------#
    # Estimating power under single kernel test
    out.row.col1 <- Single.MMD(m.x,n.y,resamp,kernel.choice[1], n.iter)
    cat("Single Kernel-1 in iteration",c(out.row.col1,k),"\n")
    
    out.row.col2 <- Single.MMD(m.x,n.y,resamp,kernel.choice[2], n.iter)
    cat("Single Kernel-2 in iteration",c(out.row.col2,k),"\n")
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    # Estimating power under multiple kernel test
    out.row.col3 <- Multi.MMD(m.x,n.y,resamp,kernel.choice[3], n.iter)
    cat("Multiple Kernel-1 in iteration",c(out.row.col3,k),"\n")
    
    out.row.col4 <- Multi.MMD(m.x,n.y,resamp,kernel.choice[4], n.iter)
    cat("Multiple Kernel-2 in iteration",c(out.row.col4,k),"\n")
    
    out.row.col5 <- Multi.MMD(m.x,n.y,resamp,kernel.choice[5], n.iter)
    cat("Multiple Kernel-3 in iteration",c(out.row.col5,k),"\n")
    #--------------------------------------------------------------------------#
    # Concatenating all the outputs
    out.row <- c(k,out.row.col1, out.row.col2, out.row.col3, 
                 out.row.col4, out.row.col5)
    #--------------------------------------------------------------------------#
  }
  stopCluster(cl)
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Set Choice","Single Kernel-1", "Single Kernel-2", 
                             "Multiple Kernel-1", "Multiple Kernel-2",
                             "Multiple Kernel-3")
  
  return (out.compare)
}
