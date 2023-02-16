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
d <- 150
# probability of mixture
p <- seq(0, 1, length.out = 6)
# parameter for sigma0 matrix generation
sigma.param <- 0.5
# parameter for sigma1 = c*sigma0 matrix generation
sigma.mult <- 1.25
# parameter for mean value for H1 distribution
mu.param <- 0


# Runnuing experiments
# storing power values for each dimensions
out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n, sigma.param = sigma.param, sigma.mult = sigma.mult,
                        mu.param = mu.param, d = d, p.seq = p, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

# Storing estimated power values

single.FR.d1 <- 2*(1:(n.rep))
power.FR.mat1 <- matrix(0, nrow = length(p), ncol = n.rep)

for (k in 1:length(p)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
}

power.FR.mat1 <- cbind(p,power.FR.mat1)

write.csv(power.FR.mat1, file = "Power-FR.csv")
