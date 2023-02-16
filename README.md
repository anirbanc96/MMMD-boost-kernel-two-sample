# Code for Boosting the power of kernel two sample test

This github repository contains the codes for experiments presented in our paper "Include Paper Here". We provide codes for generating Figures $1-9$ in our paper. The images and the corresponding codes can be found as follows:

* Figure $1$ in [Part 1(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Over%20sample%20size%20-%20dim2/TYPE-I%20ERROR/Results/GaussianScaleTypeISampleSize.pdf) and [Part 1(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Over%20sample%20size%20-%20dim2/POWER/Results/GaussianScalePowerSampleSize.pdf) with codes given in [Code for Type-I error](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Over%20sample%20size%20-%20dim2/TYPE-I%20ERROR/Code) and [Code for Power](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Over%20sample%20size%20-%20dim2/POWER/Code) respectively.
* Figure $2$ in [Part 2(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LocationScaleNormal/LocationScaleNormal.pdf) and [Part 2(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/Scalet10/Scalet10.pdf)
* Figure $3$ in [Part 3(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/ScaleNormalt10/ScaleNormalt10.pdf) and [Part 3(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LaplaceNormal/LaplaceNormal.pdf)

Codes for Figure $2$ and Figure $3$ is given in [Code for Dependence over Dimension](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Power%20over%20Dimensions/Code) where for generating Figures $2(a), 2(b), 3(a)$ and $3(b)$ the parameters from [Parameter 2(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LocationScaleNormal/SimDetails.txt), [Parameters 2(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/Scalet10/SimDetails.txt), [Parameters 3(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/ScaleNormalt10/SimDetails.txt) and [Parameters 3(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Power%20over%20Dimensions/Results/LaplaceNormal/SimDetails.txt) should be used.

* Figure $4$ in [Part 4(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Mixture%20Alternatives/d%20%3D%2030/ScaleNormaltMixing30.pdf) and [Part 4(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Mixture%20Alternatives/d%20%3D%20150/ScaleNormaltMixing150.pdf) with codes given in [Code for Mixture Alternative](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Mixture%20Alternatives/Code) where Figure $4(a)$ and $4(b)$ are generated by setting $d=30$ and $d=150$ respectively.
* Figure $5$ in [Part 5(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/Local%20Alternatives/Results/ScaleLocalNormal.pdf) and [Part 5(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/High%20Dimensional%20Alternatives/Results/ScaleNormalDimension.pdf) with codes available in [Code for local alternatives](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/Local%20Alternatives/Code) and [Codes for high-dimensional alternatives](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/High%20Dimensional%20Alternatives/Code) respectively.
* Figures $6$ and $7$ can be found in [Noisy MNIST Images](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Results/ErrorAdded%20Images) and [Power over noise strength](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MNIST-Additive%20Noise/Results/MNISTNormalData.pdf) respectively with codes in [Noisy MNIST Image Code](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MNIST-Additive%20Noise/Code/SetImage.R) and [Code for power over noise strength](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code) respectively.
* Figure $8$ can be found in [Part 8(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Reduced%20Contrast%20and%20AWGN/Results/Set%20Digit%20Images) and [Part 8(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MNIST-Reduced%20Contrast%20and%20AWGN/Results/MNISTContrastData.pdf) with codes given in [Code for Reduced Contrast Digit Images](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MNIST-Reduced%20Contrast%20and%20AWGN/Code/MNISTSetImages.R) and [Code for Power in Reduced Contrast MNIST](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Reduced%20Contrast%20and%20AWGN/Code) respectively.

The Graph and Kernel based test experiments can be implemented by calling `Body.R` in appropriate sub-folders with the appropriate parameters. 

* Figure $9$ can be found in [Part 9(a)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MMDAggComparison/one-dim%20perturbed%20uniform/Results/MMDAggCompareDim1.pdf), [Part 9(b)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MMDAggComparison/mixing%20probability/ScaleNormaltMixing30.pdf) and [Part 9(c)](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MMDAggComparison/localpower/NormalLocalPower.pdf) with codes available in [Code for Power in One-Dimensional Perturbed Alternatives](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MMDAggComparison/one-dim%20perturbed%20uniform), [Code for Comparison in mixing alternatives](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MMDAggComparison/mixing%20probability/Body.R) and [Code for Comparison in local alternatives](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/blob/main/MMDAggComparison/localpower/Body.R) respectively.

# Requirements

# Data
Aside from various simulation setups we conduct two experiments using real data in the paper. One uses a noise added version of the MNIST dataset, which is generated while implementing the [Kernel Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Kernel%20Based%20Tests) as well as the [Graph Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Graph%20Based%20Tests).

The other experiment is conducted using a noisy version of MNIST datset introduced in [The n-MNIST handwritten digit dataset](https://csc.lsu.edu/~saikat/n-mnist/). Here, in addition to additive Gaussian noise the contrast of the images is also reduced. Specifically, the contrast range is scaled down to half and an additive Gaussian noise is introduced
with signal-to-noise ratio of 12. This emulates background clutter along with significant change in lighting conditions.
# Author

Anirban Chatterjee

Department of Statistics and Data Science\
University of Pennsylvania, USA
