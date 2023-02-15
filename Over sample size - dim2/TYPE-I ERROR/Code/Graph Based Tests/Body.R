################################################################################
################# Driver Code for Power over Dimensions ########################
################################################################################
source("Functions.R")

start <- Sys.time()
# Number of repetitions
n.rep <- 50
# Number of data points
n.seq <- c(50, 100, 200, 300, 400, 500)
# Dimension vector of data
d <- 2
# probability of mixture
p <- 0
# parameter for cov matrix generation under H0
sigma.param <- 1
# parameter for cov matrix generation under H1
sigma.mult <- 1
# parameter for mean value for H1 distribution
mu.param <- 0

# Libraries for parallelising
library(foreach)
library(doParallel)
library(snow)

cores <- strtoi(Sys.getenv("NSLOTS"))-1
cl <- makeCluster(cores, methods = FALSE, type = "MPI")

registerDoParallel(cl)

# Running experiments
out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n.seq, sigma.param = sigma.param, sigma.mult = sigma.mult,
                        mu.param = mu.param, d = d, p = p, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

stopCluster(cl)

# storing type I error values
single.FR.d1 <- 2*(1:n.rep)
power.FR.mat1 <- matrix(0, nrow = length(n.seq), ncol = n.rep)

for (k in 1:length(n.seq)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
}

power.FR.mat1 <- cbind(n.seq,power.FR.mat1)

write.csv(power.FR.mat1, file = "Power-FR.csv")
