# Library for reading .mat files
library(rmatio)

# Calling required functions from Functions.R
source("Functions.R")

# reading the reduced contrast and awgn MNIST data. For source see README.md
mnist.list <- read.mat("mnist-with-reduced-contrast-and-awgn.mat")

# modifying the data for easy control
train.x <- mnist.list$train_x/255
train.y.mat <- mnist.list$train_y

train.label <- rep(NA,dim(train.y.mat)[1])

for(i in 1:dim(train.y.mat)[1]){
  train.label[i] <- which(train.y.mat[i,] == 1)-1
}

# choosing sets corresponding to P and Q
set.choice <- vector("list", 5)

set.choice[[1]] <- rbind(c(2,4,8,9), c(3,4,7,9))
set.choice[[2]] <- rbind(c(1,2,4,8,9), c(1,3,4,7,9))
set.choice[[3]] <- rbind(c(0,1,2,4,8,9), c(0,1,3,4,7,9))
set.choice[[4]] <- rbind(c(0,1,2,4,5,8,9), c(0,1,3,4,5,7,9))
set.choice[[5]] <- rbind(c(0,1,2,4,5,6,8,9), c(0,1,3,4,5,6,7,9))

# running experiments
start <- Sys.time()
n.rep <- 30
resamp <- 150

remove(mnist.list)

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(resamp, set.choice,
                        n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

# storing the power values
single.FR.d1 <- 2*(0:(n.rep))
power.FR.mat1 <- matrix(0, nrow = length(set.choice), ncol = n.rep)

for (k in 1:length(set.choice)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
  
}

power.FR.mat1 <- cbind(c(1:5),power.FR.mat1)

write.csv(power.FR.mat1, file = "Power-FR.csv")
