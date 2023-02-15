################################################################################
########## Driver code for estimating power under Noisy MNIST ##################
################################################################################

#------------------------------------------------------------------------------#
source("Functions-GraphBased.R")
#------------------------------------------------------------------------------#
require(keras)
#------------------------------------------------------------------------------#
# Getting and exploring MNIST data

mnist <- dataset_mnist()

train.x <- mnist$train$x
train.label <- mnist$train$y

train <- matrix(0, dim(train.x)[1], 28*28)

for (i in 1:dim(train.x)[1]){
  train[i,] <- as.vector(train.x[i,,])
}

train.x <- train/255
remove(mnist)
#------------------------------------------------------------------------------#
# Choise of noise variances
error.sigma <- seq(0,1,length.out = 6)

#------------------------------------------------------------------------------#
# Running experiments

start <- Sys.time()
n.rep <- 50
resamp <- 100

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(resamp, error.sigma,
                        n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

#------------------------------------------------------------------------------#
# Saving results


single.FR.d1 <- 2*(1:(n.rep))
power.FR.mat1 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)

for (k in 1:length(error.sigma)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
  
}
power.FR.mat1 <- cbind(error.sigma,power.FR.mat1)

write.csv(power.FR.mat1, file = "Power-FR.csv")
