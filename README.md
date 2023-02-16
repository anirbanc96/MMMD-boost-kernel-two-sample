# Code for Boosting the power of kernel two sample test

This github repository contains the codes for experiments presented in our paper "Include Paper Here". We provide codes for generating Figures $1-9$ in our paper. The images and the corresponding codes can be found as follows:

* Figure $1$ in [Part 1(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Over%20sample%20size%20-%20dim2/TYPE-I%20ERROR/Results/GaussianScaleTypeISampleSize.pdf) and [Part 1(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Over%20sample%20size%20-%20dim2/POWER/Results/GaussianScalePowerSampleSize.pdf) with codes given in [Code for Part 1(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Over%20sample%20size%20-%20dim2/TYPE-I%20ERROR/Code) and [Code for Part 1(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Over%20sample%20size%20-%20dim2/POWER/Code) respectively.
* Figure $2$ in [Part 2(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LocationScaleNormal/LocationScaleNormal.pdf) and [Part 2(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/Scalet10/Scalet10.pdf)
* Figure $3$ in [Part 3(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/ScaleNormalt10/ScaleNormalt10.pdf) and [Part 3(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LaplaceNormal/LaplaceNormal.pdf)

Codes for Figure $2$ and Figure $3$ is given in [Code for Dependence over Dimension](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Power%20over%20Dimensions/Code) where for generating Figures $2(a), 2(b), 3(a)$ and $3(b)$ the parameters from [Parameter 2(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LocationScaleNormal/SimDetails.txt), [Parameters 2(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/Scalet10/SimDetails.txt), [Parameters 3(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/ScaleNormalt10/SimDetails.txt) and [Parameters 3(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LaplaceNormal/SimDetails.txt) should be used.

* Figure $4$ in [Part 4(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Mixture%20Alternatives/d%20%3D%2030/ScaleNormaltMixing30.pdf) and [Part 4(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Mixture%20Alternatives/d%20%3D%20150/ScaleNormaltMixing150.pdf)
* Figure $5$ in [Part 5(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Local%20Alternatives/Results/ScaleLocalNormal.pdf) and [Part 5(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/High%20Dimensional%20Alternatives/Results/ScaleNormalDimension.pdf)

# Requirements

# Data
Aside from various simulation setups we conduct two experiments using real data in the paper. One uses a noise added version of the MNIST dataset, which is generated while implementing the [Kernel Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Kernel%20Based%20Tests) as well as the [Graph Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Graph%20Based%20Tests).

The other experiment is conducted using a noisy version of MNIST datset introduced in [The n-MNIST handwritten digit dataset](https://csc.lsu.edu/~saikat/n-mnist/). Here, in addition to additive Gaussian noise the contrast of the images is also reduced. Specifically, the contrast range is scaled down to half and an additive Gaussian noise is introduced
with signal-to-noise ratio of 12. This emulates background clutter along with significant change in lighting conditions.
# Author

Anirban Chatterjee

Department of Statistics and Data Science\
University of Pennsylvania, USA
