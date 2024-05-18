require(tidyverse)

PowerTime <- read_csv("PCombOverMean.csv")[,2:6]
colnames(PowerTime) <- c("Dimension", "Bonferroni", "HM", "Bonferroni and GM",
                         "MMMD")
d <- (PowerTime)%>%pull(1)

pdf(file="ScaleNormalDimension.pdf")

plot(d, PowerTime$Bonferroni, type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Dimension", 
     ylab='Power', main="Gaussian Location")
points(d, PowerTime$HM, type='b', col=2, pch=2, lwd=1.5)
points(d, PowerTime$`Bonferroni and GM`, type='b', col=3, pch=3, lwd=1.5)
points(d, PowerTime$MMMD, type='b', col=4, pch=4, lwd=1.5)

legend("bottomleft", c("Bonferroni","Harmonic Mean","Bonferroni and Geometric Mean",
                       "Gauss MMMD"), 
       bg='transparent', col=c(1,2,3,4), pch = c(1, 2, 3, 4))



dev.off()