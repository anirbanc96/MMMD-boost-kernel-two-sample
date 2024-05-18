################################################################################
########################### mcsapply function ##################################
require(parallel)

mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


################################################################################
################## Computing power for different tests #########################
################################################################################

source("GenData.R")
source("AuxFunctions.R")
source("PValueCombined.R")
source("MMMDTests.R")

power.d <- function(n, sigma.param, sigma.mult, 
                    mu.param, d.seq, p,
                    kernel.choice, n.est, n.iter = 500){
  
  ##############################################################################
  # input: n.seq <- vector of sample values
  #        sigma.param <- parameter for generating sigma matrix under H0
  #        sigma.mult <- parameter for multiplying cov matrix under H0
  #        mu.param <- parameter for generating mean vector under H0
  #        d <- dimension of data
  #        p <- probability of mixing
  #        kernel.choice <- vector of strings for kernel choices for each test
  #        n.iter <- number of iterations to be done for power estimation
  # output: a dataframe with estimated power for each test
  ##############################################################################
  
  set.seed(2024)
  
  writeLines(c(""), "log.txt")
  
  d <- d.seq
  out.compare <- c()
  
  #--------------------------------------------------------------------------#
  library(LaplacesDemon)
  library(Rfast)
  
  for (k in 1:length(d)){
    
    #--------------------------------------------------------------------------#
    Sigma0 <- diag(sigma.param, d[k], d[k])
    # cov matrix under H0
    Sigma1 <- sigma.mult*Sigma0
    #--------------------------------------------------------------------------#
    
    # mean vector under H0
    mu0 <- rep(0, d[k])
    # mean vector under H1
    mu1 <- rep(mu.param/sqrt(d[k]), d[k])
    
    gen.var <- list(mu0, mu1, Sigma0, Sigma1)
    #--------------------------------------------------------------------------#
    
    # Estimating power under single kernel test
    start <- Sys.time()
    out.row.col1 <- Single.MMD(n, d[k],gen.var, p,kernel.choice[1], n.est,
                               type = "bon", n.iter)
    cat(paste("Bonferroni Combination - Power: ",out.row.col1," Dimension: ",k,"\n"),
        file="log.txt", append=TRUE)
    end <- Sys.time()
    out.row.col1.time <- difftime(end, start, units = "secs")
    
    start <- Sys.time()
    out.row.col2 <- Single.MMD(n, d[k],gen.var, p,kernel.choice[1], n.est,
                               type = "hm", n.iter)
    cat(paste("Harmonic Mean Combination - Power: ",out.row.col2," Dimension: ",k,"\n"),
        file="log.txt", append=TRUE)
    end <- Sys.time()
    out.row.col2.time <- difftime(end, start, units = "secs")
    
    start <- Sys.time()
    out.row.col3 <- Single.MMD(n, d[k],gen.var, p,kernel.choice[1], n.est,
                               type = "bongm", n.iter)
    cat(paste("Bonferroni and Geometric Mean Combination - Power: ",out.row.col3," Dimension: ",k,"\n"),
        file="log.txt", append=TRUE)
    end <- Sys.time()
    out.row.col3.time <- difftime(end, start, units = "secs")
    
    #--------------------------------------------------------------------------#
    # Estimating power under multiple kernel test
    
    start <- Sys.time()
    out.row.col4 <- Multi.MMD(n, d[k],gen.var, p,kernel.choice[2], n.est, n.iter)
    cat(paste("MMMD Test - Power: ",out.row.col4," Dimension: ",k,"\n"),
        file="log.txt", append=TRUE)
    end <- Sys.time()
    out.row.col4.time <- difftime(end, start, units = "secs")
    
    cat(paste("\n"),
        file="log.txt", append=TRUE)
    
    #--------------------------------------------------------------------------#
    # Concatenating all the outputs
    out.row <- c(d[k],out.row.col1, out.row.col2, out.row.col3,
                 out.row.col4, out.row.col1.time, out.row.col2.time,
                 out.row.col3.time, out.row.col4.time)
    out.compare <- rbind(out.compare, out.row)
    
  }
  
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Sample Size","Bonferroni", "HM", "BonGM",
                             "Multiple Gaussian Kernel", "Time - Bonferroni",
                             "Time - HM", "Time - BonGM",
                             "Time - Multiple Gaussian Kernel")
  
  return (out.compare)
}
