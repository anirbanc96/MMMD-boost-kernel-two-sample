########################### FR Test Functions ##################################

############################ Distance Matrix ###################################

FR.dist <- function(z){
  
  ##############################################################################
  # INPUT: 
  # z <- a data matrix (in this case the data matrix of combination of X
  #              and y)
  # OUTPUT:
  # A distance matrix corresponding to the above data matrix
  ##############################################################################
  
  z.dist.mat <- as.matrix(dist(z))
  
  return(z.dist.mat)
}

################################################################################

########################### Counting Duplicates ################################

count.dups <- function(DF){
  
  ##############################################################################
  # INPUT: 
  # DF <- A data matrix (in this case the combined data)
  # OUTPUT:
  # A dataframe with duplicates removed and count of repetition as well
  # as first occurrence of each data vector
  ##############################################################################
  
  DT <- data.table(DF)
  DT[, list(N= .N, Index=.I[1L]), by = names(DT)]
}


################################################################################

################# Function for conducting the Graph Based Test #################

FR.test <- function(m,n, resamp.size, test.type,n.iter = 1000){
  
  ##############################################################################
  # INPUT : 
  # m,n <- sample size of two data sets
  # resamp.size <- number of resamples to be done from each data set
  # test.type <- type of graph based test to be done. See gtests documentation
  #              for more details.
  # n.iter <- number of iterations to be done for estimating power
  # OUTPUT: Estimated power of the Graph based test
  ##############################################################################
  
  count <- 0
  
  for (i in 1:n.iter){
    
    Z.dups <- count.dups(rbind(data.X[resamp.x[i,],], data.Y[resamp.y[i,],]))
    
    Z <- as.matrix(Z.dups[,-c("N","Index")])
    
    counts.mat <- cbind(1:dim(Z.dups)[1],Z.dups$N)  
    
    X.ID <- which(Z.dups$Index <= resamp.size)
    Y.ID <- which(Z.dups$Index > resamp.size)
    
    dist.mat <- FR.dist(Z)
    similarity.mat <- gTests::getGraph(counts.mat,dist.mat,5)
    
    FRtest <- gTests::g.tests(similarity.mat, X.ID, Y.ID, 
                              test.type = test.type)
    count <- count + (FRtest[1][[1]][2]<0.05)
  }
  return (count/n.iter)
}

############################ Power of FR Test ##################################
# Function for comparing power over dimensions

power.d <- function(resamp, set.choice, n.iter = 500){
  
  ##############################################################################
  # INPUTS:
  # resamp <- number of resamples to be done
  # set.choice <- vector of sets chosen
  # n.iter <- number of iterations for estimating power and threshold
  # OUTPUT:
  # A data frame having power of FR test for different set choices
  ##############################################################################
  
  # Libraries for parallelising
  library(foreach)
  library(doParallel)
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  
  registerDoParallel(cl)
  
  # Creating log file for keeping track of iterations
  writeLines(c(""), "log.txt")
  
  out.compare <- foreach(k=1:length(set.choice), .combine=rbind, .export = ls(envir=globalenv())) %dopar% {
    
    # Loading required libraries
    library(LaplacesDemon)
    library(Rfast)
    library(SpatialPack)
    library(data.table)
    
    # choice of set for H0 and H1 
    set.x <- set.choice[[k]][1,]; set.y <- set.choice[[k]][2,]
    
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
    
    # data under H0
    data.X <- train.x[index.x,]
    # data under H1
    data.Y <- train.x[index.y,]
    
    # size of data under H0 and H1
    m.x <- dim(data.X)[1]
    n.y <- dim(data.Y)[1]
    
    # indices under every resamples which is done "resamp" many times
    resamp.x <- matrix(NA, nrow = n.iter, ncol = resamp)
    resamp.y <- matrix(NA, nrow = n.iter, ncol = resamp)
    
    for(i in 1:n.iter){resamp.x[i,] <- sample(1:m.x, resamp, rep = T)}
    for(i in 1:n.iter){resamp.y[i,] <- sample(1:n.y, resamp, rep = T)}
    
    # Writing in the log file to keep track of progress
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",k,"\n"))
    
    # Estimating power under FR test
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
