Quad <- read.csv("QuadraticTime.csv")
Lin <- read.csv("LinearTime.csv")

Quad.M1 <- Quad$Multiple.Kernel.1
Quad.M2 <- Quad$Multiple.Kernel.2
Quad.M3 <- Quad$Multiple.Kernel.3

Lin.M1 <- Lin$Multiple.Kernel.1
Lin.M2 <- Lin$Multiple.Kernel.2
Lin.M3 <- Lin$Multiple.Kernel.3

Lin.S1 <- Lin$Single.Kernel.1
Lin.S2 <- Lin$Single.Kernel.2

d <- Quad$h


pdf(file="LinQuadGaussianScaled75.pdf")

plot(d, Quad.M2, type='b', col=1, pch=1, lwd=1.5, ylim=c(0, 1), xlab="Signal Strength", 
     ylab='Power', main="Gaussian Scale Alternatives")
points(d, Quad.M1, type='b', col=2, pch=2, lwd=1.5)
points(d, Quad.M3, type='b', col=3, pch=3, lwd=1.5)
points(d, Lin.M2, type='b', col=4, pch=4, lwd=1.5)
points(d, Lin.M1, type='b', col=5, pch=5, lwd=1.5)
points(d, Lin.M3, type='b', col=6, pch=6, lwd=1.5)
points(d, Lin.S1, type='b', col=7, pch=7, lwd=1.5)
points(d, Lin.S2, type='b', col=8, pch=8, lwd=1.5)


legend(1.3019, 0.65, c("Gauss MMMD","LAP MMMD","Mixed MMMD","Gauss MLMMD",
                        "LAP MLMMD","Mixed MLMMD","Gauss LMMD", "LAP LMMD"), 
       bg='transparent', col=c(1,2,3,4,5,6,7,8),
       pch = c(1, 2, 3, 4, 5, 6, 7, 8))



dev.off()
