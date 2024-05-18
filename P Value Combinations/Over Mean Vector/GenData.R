################################################################################
##################### Functions for Generating Samples #########################
################################################################################

####################### Generate Samples under H0 ##############################

X.gen <- function(n, dim, p, gen.var){
  
  ##############################################################################
  # input: n <- number of samples (integer)
  #        dim <- dimension of the data (integer)
  #        p <- probability of mixing (real number in [0,1])
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
  # input: n <- number of samples (integer)
  #        dim <- dimension of the data (integer)
  #        p <- probability of mixing (real number in [0,1])
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
