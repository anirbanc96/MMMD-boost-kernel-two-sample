


GEXP<-read.csv("MultiPower-GEXP.csv")
LAPMultiple<-read.csv("MultiPower-LAP.csv")
Mixed<-read.csv("MultiPower-MIXED.csv") 
FR<-read.csv("Power-FR.csv")
GAUSSSingle<-read.csv("SinglePower-GAUSS.csv")
LAPLACESingle<-read.csv("SinglePower-LAP.csv")

d=GEXP[,2]



pdf(file="ScaleNormaltMixing150.pdf")

plot(d, rev(apply(GEXP[,-c(1,2)], 1, mean)), type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1.2), xlab="Mixing Probability", 
ylab='Power', main="Gaussian and t-distribution Mixture in Dimension 150")
points(d, rev(apply(LAPMultiple[,-c(1,2)], 1, mean)), type='b', col=2, pch=2, lwd=1.5)
points(d, rev(apply(Mixed[,-c(1,2)], 1, mean)), type='b', col=3, pch=3, lwd=1.5)
points(d, rev(apply(FR[,-c(1,2)], 1, mean)), type='b', col=4, pch=4, lwd=1.5)
points(d, rev(apply(GAUSSSingle[,-c(1,2)], 1, mean)), type='b', col=5, pch=5, lwd=1.5)
points(d, rev(apply(LAPLACESingle[,-c(1,2)], 1, mean)), type='b', col=6, pch=6, lwd=1.5)


legend("topleft", c("Gauss MMMD","LAP MMMD","Mixed MMMD","FR", "Gauss MMD","LAP MMD"), 
       bg='transparent', col=c(1,2,3,4,5,6), pch = c(1, 2, 3, 4, 5, 6))



dev.off()


