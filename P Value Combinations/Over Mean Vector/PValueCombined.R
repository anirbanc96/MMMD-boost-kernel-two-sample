################################################################################
########################## Functions for Single Kernel #########################
################################################################################

####################### P-value for single kernel test #########################

single.H0.pvalue <- function(n, x, k, MMD.value, n.iter = 1000){
  
  ##############################################################################
  # input: n <- number of samples
  #        x <- data under H0
  #        k <- chosen kernel (using kernlab package)
  #        MMD.value <- the value of MMD statistic from observed data
  #        n.iter <- number of iterations done to estimate quantile
  # output: the estimated p-value of the test using multiplier bootstrap
  ##############################################################################
  
  require(psych)
  
  C <- diag(1, nrow = n, ncol = n) - (1/n)*matrix(1, nrow = n, ncol = n)
  
  k.mat <- (1/n)*(C%*%kernlab::kernelMatrix(k, x)%*%C)
  
  u.mat <- MASS::mvrnorm(n.iter, mu = rep(0,n), 
                         Sigma = diag(2, nrow = n, ncol = n))
  
  test.stat <- colSums(t(u.mat) * (k.mat %*% t(u.mat))) - 2*tr(k.mat)
  
  return (sum(test.stat > MMD.value)/n.iter)
}

################################################################################

################### Estimating Power for single kernel #########################

Single.MMD <- function(n, d, gen.var, p = 0, kernel.choice, n.est, type,
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
  #        n.est <- number of iterations for multiplier bootstrap
  #        type <- type of p-value combination
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
                     
                     kernel.list <- lapply(-2:2, 
                      function(x){kernlab::rbfdot(sigma = 2^(x) * sigma.med)})
                     
                     MMD.val <- mcsapply(1:5,
                                function(x){n*compute.MMD(X,Y,kernel.list[[x]])})
                     
                     MMD.pvalues <- mcsapply(1:5,
                            function(x){single.H0.pvalue(n, X,
                                                         kernel.list[[x]],
                                                         MMD.val[x], n.est)})
                     
                     if (type == "bon"){
                       
                       (5 * min(MMD.pvalues) <= 0.05)
                     }
                     
                     else if (type == "hm"){
                       
                       (2.214749 * log(5) * (5/sum(1/MMD.pvalues)) <= 0.05)
                       
                     }
                     
                     else{
                       
                       bon = 5 * min(MMD.pvalues)
                       gm = exp(1) * exp(mean(log(MMD.pvalues)))
                       
                       (2 * min(c(bon, gm)) <= 0.05)
                       
                     }
                   }
  return (sum(count)/n.iter)
}

################################################################################