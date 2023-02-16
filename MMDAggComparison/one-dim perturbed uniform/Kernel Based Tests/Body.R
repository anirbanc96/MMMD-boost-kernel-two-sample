#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
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
# Choice of kernel for single and multiple kernel tests
# poissble choices: 
# single - "GAUSS" or "LAP" for gaussian or laplacian
# multiple - "MINMAX", "GEXP" or "MIXED" for gaussian kernel with min-max, 
# exponential bandwidth choice, or mixture of gaussian and laplace kernel and
# "LAP" for laplace kernel with exponential bandwidth.

# First two coordinate for single, last three for multiple
kernel.choice <- c("LAP","GAUSS", "LAP","GEXP", "MIXED")
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
reticulate::source_python("Agg.py")
# storing power values for each dimensions


# Libraries for parallelising
library(foreach)
library(doParallel)
#library(snow)

#cores <- strtoi(Sys.getenv("NSLOTS"))-1
#cl <- makeCluster(cores, methods = FALSE, type = "MPI")

#registerDoParallel(cl)

cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer

registerDoParallel(cl)

out.d <- c()
for (iter in 1:n.rep){
  # storing power values for particular iteration
  out.d.iter <- power.d(n = n, d = d, p.seq = p.seq, sob_smooth = sob_smooth,
                        c_d = c_d, kernel.choice = kernel.choice,
                        theta_Seed = theta_seed,  n.iter = 500)
  out.d <- c(out.d, out.d.iter)
  print (iter)
}
out.d <- as.matrix(as.data.frame(out.d))
end <- Sys.time()
end-start

stopCluster(cl)

single.power.d1 <- 2+6*(0:(n.rep-1))
single.power.d2 <- 3+6*(0:(n.rep-1))
multi.power.d1 <- 4+6*(0:(n.rep-1))
multi.power.d2 <- 5+6*(0:(n.rep-1))
multi.power.d3 <- 6*(1:n.rep)
power.single.mat1 <- matrix(0, nrow = length(p.seq), ncol = n.rep)
power.single.mat2 <- matrix(0, nrow = length(p.seq), ncol = n.rep)
power.multi.mat1 <- matrix(0, nrow = length(p.seq), ncol = n.rep)
power.multi.mat2 <- matrix(0, nrow = length(p.seq), ncol = n.rep)
power.multi.mat3 <- matrix(0, nrow = length(p.seq), ncol = n.rep)
for (k in 1:length(p.seq)){
  power.single.mat1[k,] <- out.d[k,single.power.d1]
  power.single.mat2[k,] <- out.d[k,single.power.d2]
  power.multi.mat1[k,] <- out.d[k,multi.power.d1]
  power.multi.mat2[k,] <- out.d[k,multi.power.d2]
  power.multi.mat3[k,] <- out.d[k,multi.power.d3]
}

power.single.mat1 <- cbind(p.seq,power.single.mat1)
power.single.mat2 <- cbind(p.seq,power.single.mat2)
power.multi.mat1 <- cbind(p.seq,power.multi.mat1)
power.multi.mat2 <- cbind(p.seq,power.multi.mat2)
power.multi.mat3 <- cbind(p.seq,power.multi.mat3)

write.csv(power.single.mat1, file = "SinglePower-LAP.csv")
write.csv(power.single.mat2, file = "SinglePower-GAUSS.csv")
write.csv(power.multi.mat1, file = "MultiPower-LAP.csv")
write.csv(power.multi.mat2, file = "MultiPower-GEXP.csv")
write.csv(power.multi.mat3, file = "MultiPower-MIXED.csv")
