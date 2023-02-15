################################################################################
########## Driver code for estimating power under Noisy MNIST ##################
################################################################################

#------------------------------------------------------------------------------#
source("Functions.R")
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
# Choice of Kernels
# Note: First two indices are choice for single kernel tests and last three are 
# choice for multiple kernel tests
kernel.choice <- c("LAP","GAUSS", "LAP", "GEXP","MIXED")
#------------------------------------------------------------------------------#
# Running experiments
start <- Sys.time()
n.rep <- 50
resamp <- 100

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(resamp, error.sigma,
                        kernel.choice = kernel.choice, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start
#------------------------------------------------------------------------------#
# Saving results

single.power.d1 <- 2 + 6*(0:(n.rep-1))
single.power.d2 <- 3 + 6*(0:(n.rep-1))
multi.power.d1 <- 4 + 6*(0:(n.rep-1))
multi.power.d2 <- 5 + 6*(0:(n.rep-1))
multi.power.d3 <- 6*(1:n.rep)
power.single.mat1 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)
power.single.mat2 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)
power.multi.mat1 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)
power.multi.mat2 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)
power.multi.mat3 <- matrix(0, nrow = length(error.sigma), ncol = n.rep)
for (k in 1:length(error.sigma)){
  power.single.mat1[k,] <- out.d[k,single.power.d1]
  power.single.mat2[k,] <- out.d[k,single.power.d2]
  power.multi.mat1[k,] <- out.d[k,multi.power.d1]
  power.multi.mat2[k,] <- out.d[k,multi.power.d2]
  power.multi.mat3[k,] <- out.d[k,multi.power.d3]
}

power.single.mat1 <- cbind(error.sigma,power.single.mat1)
power.single.mat2 <- cbind(error.sigma,power.single.mat2)
power.multi.mat1 <- cbind(error.sigma,power.multi.mat1)
power.multi.mat2 <- cbind(error.sigma,power.multi.mat2)
power.multi.mat3 <- cbind(error.sigma,power.multi.mat3)

write.csv(power.single.mat1, file = "SinglePower-LAP.csv")
write.csv(power.single.mat2, file = "SinglePower-GAUSS.csv")
write.csv(power.multi.mat1, file = "MultiPower-LAP.csv")
write.csv(power.multi.mat2, file = "MultiPower-GEXP.csv")
write.csv(power.multi.mat3, file = "MultiPower-MIXED.csv")
#------------------------------------------------------------------------------#