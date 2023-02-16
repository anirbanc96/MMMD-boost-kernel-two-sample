# Code for Boosting the power of kernel two sample test

This github repository contains the codes for experiments presented in our paper "Include Paper Here". We provide codes for generating Figures 1-9 in our paper. The images can be found in result sub-folder of each experiment folder. 

# Requirements

# Data
Aside from various simulation setups we conduct two experiments using real data in the paper. One uses a noise added version of the MNIST dataset, which is generated while implementing the [Kernel Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Kernel%20Based%20Tests) as well as the [Graph Based Tests](https://github.com/anirbanc96/MMMD-boost-kernel-two-sample/tree/main/MNIST-Additive%20Noise/Code/Graph%20Based%20Tests).

The other experiment is conducted using a noisy version of MNIST datset introduced in [The n-MNIST handwritten digit dataset](https://csc.lsu.edu/~saikat/n-mnist/). Here, in addition to additive Gaussian noise the contrast of the images is also reduced. Specifically, the contrast range is scaled down to half and an additive Gaussian noise is introduced
with signal-to-noise ratio of 12. This emulates background clutter along with significant change in lighting conditions.
# Author

Anirban Chatterjee
Department of Statistics and Data Science
University of Pennsylvania, USA
