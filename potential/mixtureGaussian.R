### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation of gradient and density of Gaussian mixture posterior of arbitrary number of Gaussians. 
### And two additional functions to support the calculation. 

library(Matrix)

###============= mixgaussian_par: supportive function 1==========
## calculation on the parameters of the multivariate normal, like the matrix inverse, so avoid repeating such calculation in sampling. 

## input: 
#   weight: (p1, p2, p3,..., pm), a m*1 vector of weights sums up to 1, each entry a weight. m = number of Gaussian in the model. 
#   Sigma:  array of dimension d*d*m. Sigma[,,i] is the variance matrix of the i-th Gaussian. 

## output: 
#   invsigma1, invsigma2: inverse of the covariance matrix. 
#   scale1, scale2: the scale in front of the exp of normal. 
#=======================
mixgaussian_par <- function(weight, Sigma)
{
  dimension = dim(Sigma)[1]
  # dimension of space. 
  n_mixture = dim(Sigma)[3]
  # number of Gaussians in the mixture model. 
  detSigma = apply(Sigma, n_mixture, det)
  # a vector of length m. Entry i = the determinant of the variance matrix of i-th Gaussian. 
  
  mixture_scale = weight/sqrt((2*pi)^dimension*detSigma) 
  # so that scale * exp(-matrix norm of (x - mu)/2) = a vector of weighted density of each Gaussian.  
  invSigma = apply(Sigma, n_mixture, solve)
  invSigma = array(invSigma,dim = c(dimension,dimension,n_mixture))
  # inverse of the variance matrices. 

  return(list("mixture_scale" = mixture_scale, "invSigma" = invSigma))
}



###============= calculate_density_grad: supportive function 2==========

## calculate_density_grad: a function calculate the posterior density and gradient of log density wrt a given x
## input:
#    ParResult:     The outputs of twogaussian_par. 
#    mu:  Matrix of size (d, m). Column i = mean vector of i-th Gaussian. 
#    x: The current position of a sample.  

## output:
#    pd: The posterior density at x. 
#    gradU: The gradient of -log posterior. 

#======================
calculate_density_grad <- function(ParResult, mu, x)
{
  mixture_scale = ParResult$mixture_scale 
  # a vector of length m = #Gaussians. 
  n_mixture = length(mixture_scale) 
  # number of Gaussians in the mixture model. 
  dimension = dim(mu)[1]
  invSigma = ParResult$invSigma 
  # inverse of the variance matrices. 
  
  xtomu = x - mu  
  # matrix of dimension * #Gaussian
  left_xtomu = array(xtomu, dim = c(length(xtomu),1)) 
  Sigma_wrap = bdiag(apply(invSigma, n_mixture, bdiag))
  right_xtomu = bdiag(split(xtomu, rep(1:n_mixture, each = dimension)))
  
  Gaussian_dist = exp(-t(left_xtomu) %*% Sigma_wrap %*% right_xtomu)
  # the values of the Gaussian kernel within each of the Gaussians. 
  
  weighted_pd = mixture_scale * Gaussian_dist 
  weighted_pd = as.numeric(weighted_pd)
  # weighted probability densities of individual Gaussians in the mixture model. 
  pd = sum(weighted_pd)
  # probability density at x. 

  # gradU = ( weighted_pd1 * invsigma1 %*% xtomu1 + weighted_pd2 * invsigma2 %*% xtomu2)/pd
  gradU =  array(invSigma,dim = c(2,2*3)) %*% (rep(weighted_pd,each = 2) * left_xtomu)/pd
  # gradient of U at x. U = -log(probability density)
  
  return(list("pd" = pd, "gradU" = gradU))
}


### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
## parameter setting. 
dimension = 2

weight = c(1/3,1/3,1/3)
Sigma = array(NA, dim = c(2,2,3))
mu = matrix(NA, nrow = 2, ncol = 3)
for (i in 1:3) {
  Sigma[,,i] = diag(i,2)
  mu[,i] = c(i,i)
}

ParResult = mixgaussian_par(weight, Sigma)

# ===================================
## calculate the probability density and gradient for a position x.
## input: 
#    x:   positions. matrix of p*n. p = dimension of space, n = number of chains. 
## output: 
#    a data frame of two columns
#         $grad = gradient of the potential at x. 
#         $pd   = (posterior) probability density at x. 

## side note: potential = -log(probability density)
mixgaussian <- function(x) 
{
  temp <- function(x)
  {
    temp = calculate_density_grad(ParResult, mu, x)
  }
  calc_density_grad = apply(x, 2, temp)
  grad = c()
  pd   = c()
  for(i in c(1:dim(x)[2]))
  {
    grad = cbind(grad, calc_density_grad[i][[1]]$gradU)  # gradient matrix. each column a full gradient for each chain. 
    pd   = cbind(pd, calc_density_grad[i][[1]]$pd)     # vector of probability density
  }
  result = list("grad" = grad,"pd" = pd)
  return(result)
}




# #====================== test block ==================
# weight = c(1/3, 1/3, 1/3)
# Sigma = array(NA, dim = c(2,2,3))
# mu = matrix(NA, nrow = 2, ncol = 3)
# for (i in 1:3) {
#   Sigma[,,i] = diag(i,2)
#   mu[,i] = c(i,i)
# }

# x = t(mvrnorm(10, mu = rep(0,2), Sigma = diag(1,2)))
# 
# ParResult = mixgaussian_par(weight, Sigma)
# ParResult$mixture_scale
# ParResult$invSigma
# 
# # x = matrix(c(1,1,0,0, 0,0,1,1),nrow = 4)
# result = mixgaussian(x)
# result$grad



