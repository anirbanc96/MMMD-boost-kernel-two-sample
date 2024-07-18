################################################################################
################# Driver Code for Power over Dimensions ########################
################################################################################
source("Functions.R")
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
start <- Sys.time()
# Number of repetitions
n.rep <- 1
# Number of data points
n <- 500
# Dimension vector of data
d <- 1
# number of perturbations
p.seq <- c(1,2,3,4,5,6)

theta_seed <- 1
sob_smooth <- 1
c_d <- 2.7
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# storing power values for each dimensions



out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n = n, d = d, p.seq = p.seq, sob_smooth = sob_smooth,
                        c_d = c_d, theta_Seed = theta_seed, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start


single.FR.d1 <- 2*(1:n.rep)
power.FR.mat1 <- matrix(0, nrow = length(p.seq), ncol = n.rep)

for (k in 1:length(p.seq)){
  power.FR.mat1[k,] <- out.d[k,single.FR.d1]
}

power.FR.mat1 <- cbind(p.seq,power.FR.mat1)

write.csv(power.FR.mat1, file = "Power-FR-56.csv")
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#