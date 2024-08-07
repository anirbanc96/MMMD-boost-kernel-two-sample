# Function for calling pyhton code from MMDAgg Paper.

f <- function(x,y, kernel_choice, perm) {
  e <- new.env()
  options("reticulate.engine.environment" = e)
  
  tmp_var_x <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  assign(tmp_var_x, x, envir = e)
  
  tmp_var_y <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  assign(tmp_var_y, y, envir = e)
  
  tmp_var_kernel <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  assign(tmp_var_kernel, kernel_choice, envir = e)
  
  tmp_var_perm <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  assign(tmp_var_perm, perm, envir = e)
  
  reticulate::source_python("Test-Run.py")
  res <- reticulate::py_run_string(glue::glue("py_count = mmd_split_test(100,
                                              r.{tmp_var_x}, r.{tmp_var_y},
                                              0.05, r.{tmp_var_kernel},
                                              r.{tmp_var_perm},500)"))
  options("reticulate.engine.environment" = NULL)
  return (res)  
  
}  

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

library(reticulate)

# n <- 100
n <- 200
n.iter <- 500
# Dimension vector of data
# Note: Change the following value for different dimensions
d <- 30
# probability of mixture
p <- seq(0,1,length.out = 6)
# parameter for sigma0 matrix generation
sigma.param <- 0.5
# parameter for sigma1 = c*sigma0 matrix generation
sigma.mult <- 1.25
# parameter for mean value for H1 distribution
mu.param <- 0


kernel.choice <- c("gaussian", "laplace")

library(foreach)
library(doParallel)

cores <- detectCores()
cl <- makeCluster(cores[1]-1)

registerDoParallel(cl)

writeLines(c(""), "log.txt")

start <- Sys.time()

power_out <- foreach(k=1:length(p), .combine=cbind,
                     .export = ls(envir=globalenv())) %dopar% {
  
  library(LaplacesDemon)
  library(Rfast)
  
  # Sigma0 matrix <- cov matrix under H0
  Sigma0 <- matrix(0, nrow = d, ncol = d)
  for(i in 1:d){
    for(j in 1:d){
      Sigma0[i,j] = sigma.param^(abs(i-j))
    }
  }
  # cov matrix under H0
  Sigma1 <- sigma.mult*Sigma0
  #--------------------------------------------------------------------------#
  #--------------------------------------------------------------------------#
  # mean vector under H0
  mu0 <- rep(0, d)
  # mean vector under H1
  mu1 <- rep(mu.param, d)
  
  count <- rep(0, n.iter)
  
  sink("log.txt", append=TRUE)
  
  
  for (j in 1:n.iter) {
    
    cat(paste("Starting iteration",k,"\t", j, "\n"))
    
    # generating samples from H0 and H1
    x <- X.gen(n, d, p[k])
    y <- Y.gen(n, d, p[k])
    
    kernel_choice <- kernel.choice[1]
    
    perm <- "wild bootstrap"

    count[j] <- f(x,y, kernel_choice, perm)$py_count
  }
  cat(paste(mean(count)))
  mean(count)
}
stopCluster(cl)

end <- Sys.time()
end-start

write.csv(power_out, "PowerGaussMixProbOracle.csv")

