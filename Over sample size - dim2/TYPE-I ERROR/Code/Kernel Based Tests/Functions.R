################################################################################
##################### Functions for Generating Samples #########################
################################################################################

####################### Generate Samples under H0 ##############################

X.gen <- function(n, dim, p, gen.var){
  
  ##############################################################################
  # input: n <- number of samples
  #        dim <- dimension of the data
  #        p <- probability of mixing
  #        gen.var <- list containing: mu0 <- mean under null
  #                                    mu1 <- mean under alternative
  #                                    Sigma0 <- Covariance under null
  #                                    Sigma1 <- Covariance under alternative
  # output: X.samp <- n many samples from the null distribution
  ##############################################################################
  
  mu0 <- gen.var[[1]]; mu1 <- gen.var[[2]]
  Sigma0 <- gen.var[[3]]; Sigma1 <- gen.var[[4]]
  if (p == 0){
    X.samp <- MASS::mvrnorm(n, mu = mu0, Sigma = Sigma0)
  }
  else if (p == 1){
    X.samp <- LaplacesDemon::rmvt(n, mu = mu0, S = Sigma0, df = 10)
    ############################################################################
    # Uncomment following line and comment above line for generating from
    # Laplace distribution
    ############################################################################
    #X.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
  }
  else{
    choice.vec <- rbinom(n, 1, p)
    
    gauss.samp <- MASS::mvrnorm(n, mu = mu0, Sigma = Sigma0)
    other.samp <- LaplacesDemon::rmvt(n, mu = mu0, S = Sigma0, df = 10)
    ############################################################################
    # Uncomment following line and comment above line for generating from
    # Laplace distribution
    ############################################################################
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
    
    X.samp <- (1-choice.vec) * gauss.samp + choice.vec * other.samp
  }
  return (X.samp)
}

################################################################################
###################### Generate Samples under H1 ###############################

Y.gen <- function(n, dim, p, gen.var){
  
  ##############################################################################
  # input: n <- number of samples
  #        dim <- dimension of the data
  #        p <- probability of mixing
  #        gen.var <- list containing: mu0 <- mean under null
  #                                    mu1 <- mean under alternative
  #                                    Sigma0 <- Covariance under null
  #                                    Sigma1 <- Covariance under alternative
  # output: Y.samp <- n many samples from the alternative distribution
  ##############################################################################
  
  mu0 <- gen.var[[1]]; mu1 <- gen.var[[2]]
  Sigma0 <- gen.var[[3]]; Sigma1 <- gen.var[[4]]
  if (p == 0){
    Y.samp <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma1)
  }
  else if (p == 1){
    Y.samp <- LaplacesDemon::rmvt(n, mu = mu1, S = Sigma1, df = 10)
    ############################################################################
    # Uncomment following line and comment above line for generating from
    # Laplace distribution
    ############################################################################
    #Y.samp  <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
  }
  else{
    choice.vec <- rbinom(n, 1, p)
    
    gauss.samp <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma1)
    other.samp <- LaplacesDemon::rmvt(n, mu = mu1, S = Sigma1, df = 10)
    ############################################################################
    # Uncomment following line and comment above line for generating from
    # Laplace distribution
    ############################################################################
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
    
    Y.samp <- (1-choice.vec) * gauss.samp + choice.vec * other.samp
  }
  return (Y.samp)
}

################################################################################
#################### Functions for computing MMD values ########################

########################### Compute MMD_{u}^{2} ################################

compute.MMD <- function(X, Y, k){
  
  ##############################################################################
  # input: X,Y <- given data matrix
  #        k <- considerd kernel (kernlab object)
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
  #        k <- considerd kernel list (list of kernlab objects)
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
  #        k.vec <- list of kernels considered (list of kernlab objects)
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

################################################################################
########################## Functions for Single Kernel #########################
################################################################################

##################### Cutoff for single kernel under H0 ########################

single.H0.cutoff <- function(n, x, k, n.iter = 1000){
  
  ##############################################################################
  # input: n <- number of samples
  #        x <- data under H0
  #        k <- chosen kernel (kernlab object)
  #        n.iter <- number of iterations done to estimate quantile
  # output: H0.thresh <- upper alpha level threshold
  ##############################################################################
  
  require(psych)
  
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  k.mat <- (1/n)*(C%*%kernlab::kernelMatrix(k, x)%*%C)
  
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}

################################################################################

################### Estimating Power for single kernel #########################

Single.MMD <- function(n, d, gen.var, p = 0, kernel.choice = "GAUSS",
                       n.iter = 1000){
  
  ##############################################################################
  # input: n <- number of samples
  #        d <- dimension of the samples
  #        gen.var <- list containing: mu0 <- mean under null
  #                                    mu1 <- mean under alternative
  #                                    Sigma0 <- Covariance under null
  #                                    Sigma1 <- Covariance under alternative
  #        p <- mixing probability
  #        kernel.choice <- a string for choice of kernels ("GAUSS" or "LAP")
  #        n.iter <- number of iterations done for estimating power
  # output: proportion of rejections
  ##############################################################################
  
  mu0 <- gen.var[[1]]; mu1 <- gen.var[[2]]
  Sigma0 <- gen.var[[3]]; Sigma1 <- gen.var[[4]]
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p, gen.var)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, gen.var)})
  
  count <- foreach(k=1:n.iter, .combine= 'c',
                   .export = ls(envir=globalenv())) %dopar% {
    
    X <- X.list[[k]]
    Y <- Y.list[[k]]
    
    sigma.med <- med.bandwidth(X, Y)
    
    if (kernel.choice == "GAUSS"){
      kernel <- kernlab::rbfdot(sigma = sigma.med)
    }
    else if (kernel.choice == "LAP"){
      kernel <- kernlab::laplacedot(sigma = sqrt(sigma.med))
    }
    
    MMD.threshold <- single.H0.cutoff(n, X, kernel, n.iter)
    
    MMD.val <- n*compute.MMD(X,Y,kernel)
    
    (MMD.val > MMD.threshold)
    
    
  }
  return (sum(count)/n.iter)
}
################################################################################

################################################################################
########################## Functions for Mutli Kernel ##########################
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

expo.band <- function(X,Y, l0, l1){
  
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
##################### Function for multi-kernel test stat ######################

multi.func <- function(x, param){
  
  ##############################################################################
  # input: x <- a vector
  #        param <- matrix of appropriate dimensions
  # output: the value of test statistic
  ##############################################################################
  
  out <- t(x)%*%param%*%x
  
  return (out)
}

################################################################################

##################### Cutoff for multi kernel under H0 #########################

multi.k.approx.stat <- function(k.mat, u.mat){
  
  ##############################################################################
  # input: k.mat <- centered kernel matrix
  #        u.mat <- matrix of gaussian samples
  # output: value of coordinate of approximating statistic value corresponding
  #         to the given kernel
  ##############################################################################
  
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  return (test.stat)
}

################################################################################

multi.H0.cutoff <- function(n, x, k.vec, n.iter = 1000){
  
  ##############################################################################
  # input: n <- number of samples
  #        x <- dataset under H0
  #        k.vec <- list of considered kernels
  #        n.iter <- number of iterations needed to estimate cut-off
  # output: H0.thresh <- upper alpha level threshold
  ##############################################################################
  
  require(psych)
  
  k.len <- length(k.vec)
  
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  kvec.mat <- vector("list", k.len)
  for (i in 1:k.len){
    kvec.mat[[i]] <- (1/n)*(C%*%kernlab::kernelMatrix(k.vec[[i]], x)%*%C)
  }
  invcov.mat.est <- spdinv(est.cov(n, x, k.vec))
  
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  test.kernel.mat <- sapply(kvec.mat, multi.k.approx.stat, u.mat = u.mat)
  test.stat <- apply(test.kernel.mat, 1, multi.func, param = invcov.mat.est)
  
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}
################################################################################

################## Choice of kernels for multi kernel test #####################

k.choice <- function(X,Y, kernel.choice){
  
  ##############################################################################
  # input: X,Y <- given datasets
  #        kernel.choice <- string for choice of kernels 
  #                         ("MINMAX", "GEXP", "MIXED", "LAP")
  # output: list of kernels according to given choice (a list of kernlab object)
  ##############################################################################
  
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
  
  return (kernel.vec)
}

################################################################################

################### Estimating Power for multi kernel ##########################

Multi.MMD <- function(n, d,gen.var, p, kernel.choice,
                      n.iter = 1000){
  
  ##############################################################################
  # input: n <- number of samples
  #        d <- dimension of the samples
  #        gen.var <- list containing: mu0 <- mean under null
  #                                    mu1 <- mean under alternative
  #                                    Sigma0 <- Covariance under null
  #                                    Sigma1 <- Covariance under alternative
  #        p <- mixing probability
  #        kernel.choice <- a string for choice of kernels
  #                         ("MINMAX", "GEXP", "MIXED", "LAP")
  #        n.iter <- number of iterations done for estimating power
  # output: proportion of rejections
  ##############################################################################
  
  mu0 <- gen.var[[1]]; mu1 <- gen.var[[2]]
  Sigma0 <- gen.var[[3]]; Sigma1 <- gen.var[[4]]
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p, gen.var)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, gen.var)})
  
  count <- foreach(i=1:n.iter,.inorder = T,.combine = 'c',
                   .export = ls(envir=globalenv()),
                   .packages = "Rfast") %dopar% {
    
    X <- X.list[[i]]
    Y <- Y.list[[i]]

    kernel.vec <- k.choice(X, Y, kernel.choice)
    
    MMD.func.threshold <- multi.H0.cutoff(n, X, kernel.vec, n.iter)
    MMD.samp.val <- n*compute.MMD.vec(X, Y, kernel.vec)
    inv.cov.samp <- Rfast::spdinv(est.cov(n,X,kernel.vec))
    MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
    
    
    MMD.samp.func > MMD.func.threshold
    
  }
  return (sum(count)/n.iter)
}
################################################################################

################################################################################
############### Computing Type I error for different tests #####################
################################################################################

power.d <- function(n.seq, sigma.param, sigma.mult, mu.param, d, p,
                    kernel.choice, n.iter = 500){
  
  ##############################################################################
  # input: n.seq <- vector of sample values
  #        sigma.param <- parameter for generating sigma matrix under H0
  #        sigma.mult <- parameter for multiplying cov matrix under H0
  #        mu.param <- parameter for generating mean vector under H0
  #        d <- dimension of data
  #        p <- probability of mixing
  #        kernel.choice <- vector of strings for kernel choices for each test
  #        n.iter <- number of iterations to be done for type I error estimation
  # output: a dataframe with estimated type I error for each test
  ##############################################################################
  
  
  writeLines(c(""), "log.txt")
  
  n <- n.seq
  out.compare <- c()
  for (k in 1:length(n)){
    
    # Loading required libraries
    library(LaplacesDemon)
    library(Rfast)
    
    # cov matrix under H0
    Sigma0 <- diag(sigma.param, d, d)
    # cov matrix under H1
    Sigma1 <- sigma.mult*Sigma0
    
    # mean vector under H0
    mu0 <- rep(0, d)
    # mean vector under H1
    mu1 <- rep(mu.param, d)
    
    gen.var <- list(mu0, mu1, Sigma0, Sigma1)
    
    # Estimating power under single kernel test
    
    out.row.col1 <- Single.MMD(n[k], d,gen.var, p,kernel.choice[1], n.iter)
    cat(paste("Single Kernel-1 in iteration",out.row.col1," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    out.row.col2 <- Single.MMD(n[k], d,gen.var, p,kernel.choice[2], n.iter)
    cat(paste("Single Kernel-2 in iteration",out.row.col2," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    # Estimating power under multiple kernel test
    
    out.row.col3 <- Multi.MMD(n[k], d,gen.var, p,kernel.choice[3], n.iter)
    cat(paste("Multiple Kernel-1 in iteration",out.row.col3," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    out.row.col4 <- Multi.MMD(n[k], d,gen.var, p,kernel.choice[4], n.iter)
    cat(paste("Multiple Kernel-2 in iteration",out.row.col4," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    out.row.col5 <- Multi.MMD(n[k], d,gen.var, p,kernel.choice[5], n.iter)
    cat(paste("Multiple Kernel-3 in iteration",out.row.col5," ",k,"\n"),
        file="log.txt", append=TRUE)
    
    # Concatenating all the outputs
    out.row <- c(n[k],out.row.col1, out.row.col2, out.row.col3, out.row.col4,
                 out.row.col5)
    
    out.compare <- rbind(out.compare, out.row)
    
  }
  
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Sample Size","Single Kernel-1","Single Kernel-2",
                             "Multiple Kernel-1", "Multiple Kernel-2",
                             "Multiple Kernel-3")
  
  return (out.compare)
}
