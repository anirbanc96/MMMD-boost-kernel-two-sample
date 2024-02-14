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
n <- 200
# Dimension vector of data
d <- 300
# probability of mixture
#------------------------------------------------------------------------------#
# Note: Change the following values according to the simulation
p <- 0
# parameter for sigma0 matrix generation
sigma.param <- 1
# parameter for sigma1 = c*sigma0 matrix generation
sigma.mult <- seq(1,1.2,length.out = 11)
# parameter for mean value for H1 distribution
mu.param <- 0
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Choice of kernel for single and multiple kernel tests
# poissble choices: 
# single - "GAUSS" or "LAP" for gaussian or laplacian
# multiple - "MINMAX", "GEXP" or "MIXED" for gaussian kernel with min-max, 
# exponential bandwidth choice, or mixture of gaussian and laplace kernel.

# 1st three coordinate for multiple, Last two for single
kernel.choice <- c("LAP","GEXP", "MIXED", "LAP", "GAUSS")
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# storing power values for each dimensions

# Libraries for parallelising
library(foreach)
library(doParallel)
# library(snow)
# 
# cores <- strtoi(Sys.getenv("NSLOTS"))-1
# cl <- makeCluster(cores, methods = FALSE, type = "MPI")
# 
# registerDoParallel(cl)

cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer

registerDoParallel(cl)

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n, sigma.param = sigma.param, sigma.mult = sigma.mult,
                        mu.param = mu.param, d = d, p = p,
                        kernel.choice = kernel.choice, n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}

# Making the output a matrix for ease of access
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

stopCluster(cl)

write.csv(out.d, file = "Lineard300mult1.2.csv")