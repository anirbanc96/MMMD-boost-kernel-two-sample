################################################################################
########################## Functions for Mutli Kernel ##########################
################################################################################


##################### Function for multi-kernel test stat ######################

multi.func <- function(x, param){
  
  ##############################################################################
  # input: x <- dataset under H0
  #        param <- the inverted covariance matrix
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
  
  test.kernel.mat <- mcsapply(kvec.mat, multi.k.approx.stat, u.mat = u.mat)
  test.stat <- apply(test.kernel.mat, 1, multi.func, param = invcov.mat.est)
  
  H0.thresh <- quantile(test.stat, probs = 0.95)
  return (H0.thresh)
}

################################################################################


################### Estimating Power for multi kernel ##########################

Multi.MMD <- function(n, d,gen.var, p = 0, kernel.choice = "GEXP", n.est,
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
                     
                     MMD.func.threshold <- multi.H0.cutoff(n, X, kernel.vec, n.est)
                     
                     MMD.samp.val <- n*compute.MMD.vec(X, Y, kernel.vec)
                     inv.cov.samp <- Rfast::spdinv(est.cov(n,X,kernel.vec))
                     MMD.samp.func <- multi.func(MMD.samp.val, param = inv.cov.samp)
                     
                     
                     MMD.samp.func > MMD.func.threshold
                     
                   }
  return (sum(count)/n.iter)
}

################################################################################
