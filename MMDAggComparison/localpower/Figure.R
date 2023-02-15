


GEXP<-read.csv("MultiPower-GEXP.csv")
LAPMultiple<-read.csv("MultiPower-LAP.csv")
Mixed<-read.csv("MultiPower-MIXED.csv") 
FR<-read.csv("Power-FR.csv")

d=GEXP[,2]-1


library(tidyverse)
MMDAgg.gaussian.power <- read_csv("PowerGaussLocal.csv") %>%
  as.matrix()

MMDAgg.laplace.power <- read_csv("PowerLapLocal.csv") %>%
  as.matrix()



pdf(file="NormalLocalPower.pdf")

plot(d, apply(GEXP[,-c(1,2)], 1, mean), type='b', col=1, pch=1, lwd=1.5,
     ylim=c(0, 0.95), xlab="Signal Strength", 
     ylab='Power', main="Gaussian Scale Local Alternatives")
points(d, apply(LAPMultiple[,-c(1,2)], 1, mean), type='b', col=2, pch=2, lwd=1.5)
points(d, apply(Mixed[,-c(1,2)], 1, mean), type='b', col=3, pch=3, lwd=1.5)
points(d, apply(FR[,-c(1,2)], 1, mean), type='b', col=4, pch=4, lwd=1.5)
points(d, MMDAgg.gaussian.power[,-1], type='b', col=5, pch=5, lwd=1.5)
points(d, MMDAgg.laplace.power[,-1], type='b', col=6, pch=6, lwd=1.5)


legend("topleft", c("Gauss MMMD","LAP MMMD","Mixed MMMD","FR", "Gauss MMDAgg","LAP MMDAgg"), 
       bg='transparent', col=c(1,2,3,4,5,6), pch = c(1, 2, 3, 4, 5, 6))



dev.off()

