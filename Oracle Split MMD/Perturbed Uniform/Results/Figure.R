GEXP<-read.csv("MultiPower-GEXP.csv")
LAPMultiple<-read.csv("MultiPower-LAP.csv")
Mixed<-read.csv("MultiPower-MIXED.csv") 
# FR<-read.csv("Power-FR.csv")
#GAUSSSingle<-read.csv("SinglePower-GAUSS.csv")
#LAPLACESingle<-read.csv("SinglePower-LAP.csv")

d=GEXP[,2]

library(tidyverse)
results <- read_csv("results.csv")
split.gauss.power <- results %>%
  filter(experiment == "3: uniform alternative") %>%
  filter(kernel_type == "gaussian") %>%
  filter(function_type == "split") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()

oracle.gauss.power <- results %>%
  filter(experiment == "3: uniform alternative") %>%
  filter(kernel_type == "gaussian") %>%
  filter(function_type == "split (doubled sample sizes)") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()

split.lap.power <- results %>%
  filter(experiment == "3: uniform alternative") %>%
  filter(kernel_type == "laplace") %>%
  filter(function_type == "split") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()

oracle.lap.power <- results %>%
  filter(experiment == "3: uniform alternative") %>%
  filter(kernel_type == "laplace") %>%
  filter(function_type == "split (doubled sample sizes)") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()

pdf(file="SplitOracleCompareDim1.pdf")

plot(d, GEXP[,3], type='b', col=1, pch=1, lwd=1.5,
     ylim=c(0, 1), xlab="Perturbations", 
     ylab='Power', main="Perturbed One Dimensional Uniform Distribution")
points(d, LAPMultiple[,3], type='b', col=2, pch=2, lwd=1.5)
points(d, Mixed[,3], type='b', col=3, pch=3, lwd=1.5)
# points(d, FR[,3], type='b', col=4, pch=4, lwd=1.5)
points(d, split.gauss.power[,2], type='b', col=4, pch=4, lwd=1.5)
points(d, split.lap.power[,2], type='b', col=5, pch=5, lwd=1.5)
points(d, oracle.gauss.power[,2], type='b', col=6, pch=6, lwd=1.5)
points(d, oracle.lap.power[,2], type='b', col=7, pch=7, lwd=1.5)


legend("topright", c("Gauss MMMD","Lap MMMD","Mixed MMMD","Opt Gauss MMD",
                     "Opt Lap MMD", "Oracle Gauss MMD", "Oracle Lap MMD"), 
       bg='transparent', col=c(1,2,3,4,5,6,7),
       pch = c(1, 2, 3, 4, 5, 6, 7))



dev.off()


