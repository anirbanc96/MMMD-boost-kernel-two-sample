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
# Function for generating samples under alternative

Y.gen <- function(n, dim, num_perturb, sob_smooth = 1, c_d = 2.7, 
                  theta_seed = 100){
  theta_seed <- sample(1:(1000*n), 1)
  samp_seed <- sample(1:(1000*n), 1)
  Y.samp <- f_theta_sampler(theta_seed, samp_seed, n, num_perturb, sob_smooth,
                            c_d, dim)
  
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
# n.iter <- number of iterations to be done
# OUTPUTS: Estimated power of Graph Based test

FR.test <- function(n,d,p, sob_smooth = 1, c_d = 2.7, 
                    theta_seed = 1, test.type,n.iter = 1000){
  
  count <- 0
  
  counts.mat <- cbind(1:(2*n), rep(1,2*n))
  
  X.list <- lapply(1:n.iter,function(x){X.gen(n,d)})
  Y.list <- lapply(1:n.iter,function(x){Y.gen(n,d,p, sob_smooth, c_d, 
                                              theta_seed)})
  
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
################################################################################

################################################################################
########################## Power comparison function ###########################
################################################################################
# Function for comparing power over dimensions
# INPUTS:
# n <- no of data points from both distributions
# sigma.param <- the parameter used for generating variance covariance matrix
# sigma.mult <- the parameter for multiplying cov matrix under H0 to get cov
# under H1
# mu.param <- the parameter for generating mean under H1
# d.seq <- the vector of dimensions
# p <- the mixing probability
# OUTPUT:
# A data frame having powers of graph based test under various dimensions
power.d <- function(n,d, p.seq,sob_smooth,c_d,theta_Seed = 1,  n.iter = 500){
  
  # Libraries for parallelising
  library(foreach)
  library(doParallel)
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  
  # redefining dimension vector for ease
  p <- p.seq
  
  out.compare <- foreach(k=1:length(p), .combine=rbind, .export = ls(envir=globalenv())) %dopar% {
    
    #--------------------------------------------------------------------------#
    library(Rfast)
    reticulate::source_python("Agg.py")
    #--------------------------------------------------------------------------#
    # Creating a log file to keep track of progress
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",k,"\n"))
    #--------------------------------------------------------------------------#
    
    #----------------------------------------------------------------------------#
    #----------------------------------------------------------------------------#
    
    # Estimating asymptotic power under FR test
    out.row.col1 <- FR.test(n, d, p[k],sob_smooth, c_d, 
                            theta_seed, test.type = "o", n.iter)
    cat("FR test in iteration",c(out.row.col1,k),"\n")
    
    # Concatenating all the outputs
    out.row <- c(p[k],out.row.col1)
  }
  stopCluster(cl)
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Dimension", "FR test")
  
  return (out.compare)
}