#------------------------------------------------------------------------------#
#----------------------------- Graph Based Test -------------------------------#
#------------------------------------------------------------------------------#

############################ Distance Matrix ###################################\

# INPUT : z <- a data matrix (in this case the data matrix of combination of X
#              and y)

# OUTPUT : A distance matrix corresponding to the above data matrix

FR.dist <- function(z){
  
  z.dist.mat <- as.matrix(dist(z))
  
  return(z.dist.mat)
}

################################################################################

########################### Counting Duplicates ################################

# INPUT: DF <- A data matrix (in this case the combined data)
# OUTPUT <- A dataframe with duplicates removed and count of repetition as well
#           as first occurrence of each data vector

count.dups <- function(DF){
  
  DT <- data.table(DF)
  DT[, list(N= .N, Index=.I[1L]), by = names(DT)]
}

################################################################################

################# Function for conducting the Graph Based Test #################

# INPUT : 
# m,n <- sample size of two data sets
# resamp.size <- number of resamples to be done from each data set
# test.type <- type of graph based test to be done. See gtests documentation
# for more details.
# n.iter <- number of iterations to be done for estimating power

# OUTPUT: Estimated power of the Graph based test

FR.test <- function(m,n, resamp.size, test.type,n.iter = 1000){
  
  # variable for counting number of rejections
  count <- 0
  
  for (i in 1:n.iter){
    # dataframe containing unique data points with count and location of 
    # repetitions
    Z.dups <- count.dups(rbind(data.X[resamp.x[i,],], data.Y[resamp.y[i,],]))
    
    # dataframe cotaining only the unique data points
    Z <- as.matrix(Z.dups[,-c("N","Index")])
    
    # matrix with repetition counts of unqiue data points
    counts.mat <- cbind(1:dim(Z.dups)[1],Z.dups$N)  
    
    # The indices corresponding to X
    X.ID <- which(Z.dups$Index <= resamp.size)
    
    # The indices corresponding to Y
    Y.ID <- which(Z.dups$Index > resamp.size)
    
    # distance matrix between the unique data points
    dist.mat <- FR.dist(Z)
    
    # edge matrix formed using above distance matrix
    similarity.mat <- gTests::getGraph(counts.mat,dist.mat,5)
    
    # conduting the graph based test
    FRtest <- gTests::g.tests(similarity.mat, X.ID, Y.ID, 
                              test.type = test.type)
    
    # adding rejections
    count <- count + (FRtest[1][[1]][2]<0.05)
    #print (c(i,count))
  }
  return (count/n.iter)
}


#------------------------------------------------------------------------------#
#---------------------------- Power plot function ----------------------------#
#------------------------------------------------------------------------------#
# Function for comparing power 
# INPUTS: 
# resamp <- number of resamples to consider
# error.sigma <- vector of error variances
# n.iter <- number of iterations to be done to estimate power

# OUTPUT: Estimated power of the Graph Based test for each error variance
power.d <- function(resamp, error.sigma, n.iter = 500){
  
  # Libraries for parallelising
  library(foreach)
  library(doParallel)
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  
  out.compare <- foreach(k=1:length(error.sigma), .combine=rbind, .export = ls(envir=globalenv())) %dopar% {
    
    #--------------------------------------------------------------------------#
    library(LaplacesDemon)
    library(Rfast)
    library(SpatialPack)
    library(data.table)
    
    set.x <- c(1,2,3); set.y <- c(1,2,8)
    
    index.y <- rep(F,length(train.label)); index.x <- rep(F,length(train.label))
    for(i in 1:length(train.label)){
      if (train.label[i] %in% set.y){
        index.y[i] <- T
      }
      if (train.label[i] %in% set.x){
        index.x[i] <- T
      }
    }
    
    data.X <- train.x[index.x,]
    data.Y <- train.x[index.y,]
    
    m.x <- dim(data.X)[1]
    n.y <- dim(data.Y)[1]
    
    data.X <- imnoise(data.X, type = "gaussian", mean=0, sd=error.sigma[k])
    data.Y <- imnoise(data.Y, type = "gaussian", mean=0, sd=error.sigma[k])
    
    resamp.x <- matrix(NA, nrow = n.iter, ncol = resamp)
    resamp.y <- matrix(NA, nrow = n.iter, ncol = resamp)
    
    for(i in 1:n.iter){resamp.x[i,] <- sample(1:m.x, resamp, rep = T)}
    for(i in 1:n.iter){resamp.y[i,] <- sample(1:n.y, resamp, rep = T)}
    #--------------------------------------------------------------------------#
    # Creating a log file to keep track of progress
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",k,"\n"))
    #--------------------------------------------------------------------------#
    
    #----------------------------------------------------------------------------#
    #----------------------------------------------------------------------------#
    
    out.row.col1 <- FR.test(m.x,n.y,resamp, test.type = "o", n.iter)
    cat("FR test in iteration",c(out.row.col1,k),"\n")
    
    # Concatenating all the outputs
    out.row <- c(k,out.row.col1)
  }
  stopCluster(cl)
  
  out.compare <- as.data.frame(out.compare)
  colnames(out.compare) <- c("Set Choice", "FR test")
  
  return (out.compare)
}