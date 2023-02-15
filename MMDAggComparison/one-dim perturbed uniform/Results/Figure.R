


GEXP<-read.csv("MultiPower-GEXP.csv")
LAPMultiple<-read.csv("MultiPower-LAP.csv")
Mixed<-read.csv("MultiPower-MIXED.csv") 
FR<-read.csv("Power-FR.csv")
#GAUSSSingle<-read.csv("SinglePower-GAUSS.csv")
#LAPLACESingle<-read.csv("SinglePower-LAP.csv")

d=GEXP[,2]

library(tidyverse)
MMDAgg.results <- read_csv("CopyOfresults.csv")
MMDAgg.gaussian.power <- MMDAgg.results %>%
  filter(experiment == "1: uniform alternative") %>%
  filter(kernel_type == "gaussian") %>%
  filter(l_minus_l_plus == "(-6, -2)") %>%
  filter(function_type == "increasing") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()

MMDAgg.laplace.power <- MMDAgg.results %>%
  filter(experiment == "1: uniform alternative") %>%
  filter(kernel_type == "laplace") %>%
  filter(l_minus_l_plus == "(-6, -2)") %>%
  filter(function_type == "increasing") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power)) %>%
  as.matrix()


pdf(file="MMDAggCompareDim1.pdf")

plot(d, GEXP[,3], type='b', col=1, pch=1, lwd=1.5,
     ylim=c(0, 1), xlab="Perturbations", 
     ylab='Power', main="Perturbed One Dimensional Uniform Distribution")
points(d, LAPMultiple[,3], type='b', col=2, pch=2, lwd=1.5)
points(d, Mixed[,3], type='b', col=3, pch=3, lwd=1.5)
points(d, FR[,3], type='b', col=4, pch=4, lwd=1.5)
points(d, MMDAgg.gaussian.power[,2], type='b', col=5, pch=5, lwd=1.5)
points(d, MMDAgg.laplace.power[,2], type='b', col=6, pch=6, lwd=1.5)


legend("topright", c("Gauss MMMD","LAP MMMD","Mixed MMMD","FR", "Gauss MMDAgg","LAP MMDAgg"), 
       bg='transparent', col=c(1,2,3,4,5,6), pch = c(1, 2, 3, 4, 5, 6))



dev.off()


