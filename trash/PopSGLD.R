setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

source("2gaussian.R")  # cite functions for parameter setting. 
# source("trash/sample_PMSGLD.R")
# library(MASS)

# ===================================
## parameter setting
dimension = 2
mu1 = rep(0,dimension)
mu2 = rep(5,dimension)
sigma1 = diag(1,dimension)
sigma2 = diag(1,dimension)
weight = 1/3
ParResult = twogaussian_par(weight, sigma1, sigma2)


# ===================================


## calculate the probability density and gradient for each chain. 
# twogaussian function in to uni-input, in order to try using apply. 
twogaussian_density_grad <- function(x) 
{
  twogaussian_density_grad = twogaussian(ParResult, mu1, mu2, x)
}

# ===================================
## initialization of positions. 
nchain = 10
x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
x = t(x)
# x =  matrix(c(1,1,0,0, 0,0,1,1),nrow = dimension)


# ===================================
## sample according to population momentum sgld sampler

sample_path_p = sample_PMSGLD(50000, x,m)