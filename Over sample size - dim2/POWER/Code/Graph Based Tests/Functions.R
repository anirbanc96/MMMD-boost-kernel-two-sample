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
    
### Uncomment line below and comment line above for generating from Laplace ###
    
    #X.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
  }
  else{
    
    # vector of choices for mixture distributions
    choice.vec <- rbinom(n, 1, p)
    
    # generating the samples from Gaussian distribution
    gauss.samp <- MASS::mvrnorm(n, mu = mu0, Sigma = Sigma0)
    
    # generating the samples from t (df = 10) or Laplace for mixture
    other.samp <- LaplacesDemon::rmvt(n, mu = mu0, S = Sigma0, df = 10)
    
### Uncomment line below and comment line above for generating from Laplace ###
    
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu0, Sigma = Sigma0)
    
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
    
### Uncomment line below and comment line above for generating from Laplace ###
    
    #Y.samp  <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
  }
  else{
    
    # vector of choices for mixture distributions
    choice.vec <- rbinom(n, 1, p)
    
    # generating the samples from Gaussian distribution
    gauss.samp <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma1)
    
    # generating the samples from t (df = 10) or Laplace for mixture
    other.samp <- LaplacesDemon::rmvt(n, mu = mu1, S = Sigma1, df = 10)
    
### Uncomment line below and comment line above for generating from Laplace ###
    
    #other.samp <- LaplacesDemon::rmvl(n, mu = mu1, Sigma = Sigma1)
    
    Y.samp <- (1-choice.vec) * gauss.samp + choice.vec * other.samp
  }
  return (Y.samp)
}
################################################################################


################################################################################
################# Function for Power for Graph Based tests #####################
################################################################################

############### Function for computing Euclidean Distance Matrix ###############
# INPUTS:
# x <- data under H0
# y <- data under H1
# OUTPUT:  a distance matrix from the combined data x and y
FR.dist <- function(x,y){
  z <- rbind(x,y)
  z.dist.mat <- as.matrix(dist(z))
  
  return(z.dist.mat)
}


############### Function for estimating power of Graph Based tests #############
# INPUTS:
# n <- number of data points
# d <- dimension of the data
# p <- probability of mixing
# test.type <- type of test to be done. See gtests documentation for details
# n.iter <- number of iterations to be done
# OUTPUTS: Estimated power of Graph Based test

FR.test <- function(n,d,p, test.type,n.iter = 1000){
  
  count <- 0
  
  counts.mat <- cbind(1:(2*n), rep(1,2*n))
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d,p)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p)})
  
  for (i in 1:n.iter){
    
    dist.mat <- FR.dist(X.list[[i]],Y.list[[i]])
    
    ############################ MST SIZE ######################################
    # Here we have taken 5-MST 
    ############################################################################
    similarity.mat <- gTests::getGraph(counts.mat,dist.mat,5)
    FRtest <- gTests::g.tests(similarity.mat, 1:n, (n+1):(2*n), 
                              test.type = test.type)
    count <- count + (FRtest[1][[1]][2]<0.05)
    #print (c(i,count))
  }
  return (count/n.iter)
}


################################################################################
########################## Power comparison function ###########################
################################################################################
# Function for comparing power over dimensions
# INPUTS:
# n.seq <- vector of sample sizes
# sigma.param <- the parameter used for generating cov matrix under H0
# sigma.mult <- the parameter for multiplying cov matrix under H0 to get cov
# under H1
# mu.param <- the parameter for generating mean under H1
# d <- dimension of data
# p <- the mixing probability
# n.iter <- number of iterations performed for estimating power
# OUTPUT:
# A data frame having powers of graph based test for different sample sizes
power.d <- function(n.seq, sigma.param = 0.4, sigma.mult = 1.1, 
                    mu.param = 0, d,p = 0, n.iter = 500){
  
  
  writeLines(c(""), "log.txt")
  
  # redefining sample size vector for ease
  n <- n.seq
  
  out.compare <- foreach(k=1:length(n), .combine=rbind, .export = ls(envir=globalenv())) %dopar% {
    
    #--------------------------------------------------------------------------#
    library(LaplacesDemon)
    library(Rfast)
    #--------------------------------------------------------------------------#
    # Creating a log file to keep track of progress
    cat(paste("\n","Starting iteration",k,"\n"), 
       file="log.txt", append=TRUE)
    #--------------------------------------------------------------------------#
    
    # Sigma0 matrix <- cov matrix under H0
    Sigma0 <- diag(sigma.param, d, d)
    # cov matrix under H1
    Sigma1 <- sigma.mult*Sigma0
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    # mean vector under H0
    mu0 <- rep(0, d)
    # mean vector under H1
    mu1 <- rep(mu.param, d)
    
    #----------------------------------------------------------------------------#
    #----------------------------------------------------------------------------#
    
    # Estimating asymptotic power under FR test
    out.row.col1 <- FR.test(n[k], d, p, test.type = "o", n.iter)
    cat(paste("FR in iteration",out.row.col1," ",k,"\n"), file="log.txt", append=TRUE)
    
    # Concatenating all the outputs
    out.row <- c(n[k],out.row.col1)
  }
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Sample Size", "FR test")
  
  return (out.compare)
}
