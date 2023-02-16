################################################################################
########## Driver code for estimating power under Noisy MNIST ##################
################################################################################

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
#------------------------------------------------------------------------------#
kernel.choice <- c("LAP","GAUSS", "LAP", "GEXP","MIXED")
#------------------------------------------------------------------------------#
# Running Experiments
start <- Sys.time()
n.rep <- 30
resamp <- 150

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(resamp, set.choice,
                        kernel.choice = kernel.choice, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start
#------------------------------------------------------------------------------#
# Storing Results
single.power.d1 <- 2 + 6*(0:(n.rep-1))
single.power.d2 <- 3 + 6*(0:(n.rep-1))
multi.power.d1 <- 4 + 6*(0:(n.rep-1))
multi.power.d2 <- 5 + 6*(0:(n.rep-1))
multi.power.d3 <- 6*(1:n.rep)
power.single.mat1 <- matrix(0, nrow = length(set.choice), ncol = n.rep)
power.single.mat2 <- matrix(0, nrow = length(set.choice), ncol = n.rep)
power.multi.mat1 <- matrix(0, nrow = length(set.choice), ncol = n.rep)
power.multi.mat2 <- matrix(0, nrow = length(set.choice), ncol = n.rep)
power.multi.mat3 <- matrix(0, nrow = length(set.choice), ncol = n.rep)
for (k in 1:length(set.choice)){
  power.single.mat1[k,] <- out.d[k,single.power.d1]
  power.single.mat2[k,] <- out.d[k,single.power.d2]
  power.multi.mat1[k,] <- out.d[k,multi.power.d1]
  power.multi.mat2[k,] <- out.d[k,multi.power.d2]
  power.multi.mat3[k,] <- out.d[k,multi.power.d3]
}

power.single.mat1 <- cbind(1:5,power.single.mat1)
power.single.mat2 <- cbind(1:5,power.single.mat2)
power.multi.mat1 <- cbind(1:5,power.multi.mat1)
power.multi.mat2 <- cbind(1:5,power.multi.mat2)
power.multi.mat3 <- cbind(1:5,power.multi.mat3)

write.csv(power.single.mat1, file = "SinglePower-LAP.csv")
write.csv(power.single.mat2, file = "SinglePower-GAUSS.csv")
write.csv(power.multi.mat1, file = "MultiPower-LAP.csv")
write.csv(power.multi.mat2, file = "MultiPower-GEXP.csv")
write.csv(power.multi.mat3, file = "MultiPower-MIXED.csv")
