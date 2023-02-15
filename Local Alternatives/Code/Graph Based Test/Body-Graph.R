################################################################################
################# Driver Code for Power over Dimensions ########################
################################################################################
source("Functions.R")

start <- Sys.time()
# Number of repetitions
n.rep <- 1
# Number of data points
n <- 100
# Dimension vector of data
d <- 20
# probability of mixture
p <- 0
# parameter for sigma1 = c*sigma0 matrix generation
sigma.mult <- 1+seq(0,2,length.out = 11)/sqrt(n)
# parameter for mean value for H1 distribution
mu.param <- 0

# storing power values for each signal strength
out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n, sigma.mult.seq = sigma.mult,
                        mu.param = mu.param, d = d, p = p, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start


single.FR.d1 <- 2*(1:n.rep)

power.FR.mat1 <- matrix(0, nrow = length(sigma.mult), ncol = n.rep)


for (k in 1:length(sigma.mult)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
  
}

power.FR.mat1 <- cbind(sigma.mult,power.FR.mat1)


write.csv(power.FR.mat1, file = "Power-FR.csv")
