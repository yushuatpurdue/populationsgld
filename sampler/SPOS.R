### This file is created by Yushu Huang for local SPOS experiments.
### This file runs SPOS. Algorithm see Jianyi 2020. 

library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/2gaussian.R")
distribution_sgld <- twogaussian

## ================ test =================
nchain = 10
dimension = 5
x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
x = t(x)
pd_grad = distribution_sgld(x)