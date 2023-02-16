f <- function(x,y, kernel_choice, perm, type) {
  e <- new.env()
  options("reticulate.engine.environment" = e)
  
  # create a new variable with a random name
  tmp_var_x <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  assign(tmp_var_x, x, envir = e)
  
  tmp_var_y <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  assign(tmp_var_y, y, envir = e)
  
  tmp_var_kernel <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  assign(tmp_var_kernel, kernel_choice, envir = e)
  
  tmp_var_perm <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  assign(tmp_var_perm, perm, envir = e)
  
  tmp_var_type <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  assign(tmp_var_type, type, envir = e)
  
  #tmp_var_l <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  #assign(tmp_var_l, l, envir = e)
  
  #tmp_var_u <- paste(sample(letters, 30, replace = TRUE), collapse = "")
  #message(tmp_var_name)
  #assign(tmp_var_u, u, envir = e)
  
  reticulate::source_python("Test-Run.py")
  res <- reticulate::py_run_string(glue::glue("py_count = mmdagg(100, r.{tmp_var_x}, r.{tmp_var_y}, 0.05, r.{tmp_var_kernel}, r.{tmp_var_perm}, r.{tmp_var_type}, 10,14, 500, 500, 100)"))
  options("reticulate.engine.environment" = NULL)  # unset option
  return (res)  
  
}

library(reticulate)


#------------------------------------------------------------------------------#
# Library for reading .mat files
library(rmatio)
#------------------------------------------------------------------------------#
# Noise added and reduced contrast MNIST data
mnist.list <- read.mat("mnist-with-reduced-contrast-and-awgn.mat")

# modifying the data for easy control
train.x <- mnist.list$train_x/255
train.y.mat <- mnist.list$train_y

train.label <- rep(NA,dim(train.y.mat)[1])

for(i in 1:dim(train.y.mat)[1]){
  train.label[i] <- which(train.y.mat[i,] == 1)-1
}
#------------------------------------------------------------------------------#
# choosing sets corresponding to P and Q
set.choice <- vector("list", 5)

set.choice[[1]] <- rbind(c(2,4,8,9), c(3,4,7,9))
set.choice[[2]] <- rbind(c(1,2,4,8,9), c(1,3,4,7,9))
set.choice[[3]] <- rbind(c(0,1,2,4,8,9), c(0,1,3,4,7,9))
set.choice[[4]] <- rbind(c(0,1,2,4,5,8,9), c(0,1,3,4,5,7,9))
set.choice[[5]] <- rbind(c(0,1,2,4,5,6,8,9), c(0,1,3,4,5,6,7,9))

resamp <- 150
n.iter <- 500

kernel.choice <- c("gaussian", "laplace")

library(foreach)
library(doParallel)

cores <- detectCores()
cl <- makeCluster(cores[1]-1)

registerDoParallel(cl)

writeLines(c(""), "log.txt")

power_out <- foreach(k=1:length(set.choice), .combine=cbind, .export = ls(envir=globalenv())) %dopar% {
  
  set_x <- set.choice[[k]][1, ]
  set_y <- set.choice[[k]][2, ]
  
  index_y <- rep(F, length(train.label))
  index_x <- rep(F, length(train.label))
  for (i in 1:length(train.label)) {
    if (train.label[i] %in% set_y) {
      index_y[i] <- T
    }
    if (train.label[i] %in% set_x) {
      index_x[i] <- T
    }
  }
  
  # data under H0
  data_X <- train.x[index_x, ]
  # data under H1
  data_Y <- train.x[index_y, ]
  
  # size of data under H0 and H1
  m_x <- dim(data_X)[1]
  n_y <- dim(data_Y)[1]
  
  # indices under every resamples which is done "resamp" many times
  resamp_x <- matrix(NA, nrow = n.iter, ncol = resamp)
  resamp_y <- matrix(NA, nrow = n.iter, ncol = resamp)
  
  for (i in 1:n.iter) {
    resamp_x[i, ] <- sample(1:m_x, resamp, rep = T)
  }
  for (i in 1:n.iter) {
    resamp_y[i, ] <- sample(1:n_y, resamp, rep = T)
  }
  
  count <- rep(0, n.iter)
  
  sink("log.txt", append=TRUE)
  
  
  for (j in 1:n.iter) {
    
    cat(paste("Starting iteration",k,"\t", j, "\n"))
    
    x <- data_X[resamp_x[j, ], ]
    y <- data_Y[resamp_y[j, ], ]
    
    kernel_choice <- kernel.choice[1]
    
    perm <- "wild bootstrap"
    type <- "increasing"
    count[j] <- f(x,y, kernel_choice, perm, type)$py_count
  }
  cat(paste(mean(count)))
  mean(count)
}
stopCluster(cl)
write.csv(power_out, "PowerGauss.csv")

