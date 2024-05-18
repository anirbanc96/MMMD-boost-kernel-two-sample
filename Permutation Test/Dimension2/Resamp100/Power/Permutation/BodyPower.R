################################################################################
source("Functions.R")
################################################################################

start <- Sys.time()
# Number of repetitions
n.rep <- 1
# Number of data points
n.seq <- c(200, 300, 400, 500, 600, 700)
# Dimension vector of data
d <- 2
# probability of mixture
p <- 0
# parameter for sigma0 matrix generation
sigma.param <- 1
# parameter for sigma1 = c*sigma0 matrix generation
sigma.mult <- 1.25
# parameter for mean value for H1 distribution
mu.param <- 0
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Choice of kernel for single and multiple kernel tests
# poissble choices: 
# single - "GAUSS" or "LAP" for gaussian or laplacian
# multiple - "MINMAX", "GEXP" or "MIXED" for gaussian kernel with min-max, 
# exponential bandwidth choice, or mixture of gaussian and laplace kernel.

kernel.choice <- "GEXP"
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# storing power values for each dimensions


# Libraries for parallelising
library(foreach)
library(doParallel)

# Uncomment following package if using MPI
library(snow)
 
cores <- strtoi(Sys.getenv("NSLOTS"))
cl <- makeCluster(cores, methods = FALSE, type = "MPI")
 
registerDoParallel(cl)

################################################################################
# Uncommment following line if not using MPI
# cores <- detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer

# registerDoParallel(cl)
################################################################################

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n.seq, sigma.param = sigma.param,
                        sigma.mult = sigma.mult,
                        mu.param = mu.param, d, p = p,
                        kernel.choice = kernel.choice,
                        n.perm = 100, n.est = 100, n.iter = 100)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
difftime(end,start, units = "secs")

# stopCluster(cl)
################################################################################
write.csv(out.d, file = "TimePowerCompare.csv")