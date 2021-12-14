### This file is created by Yushu Huang for population SGMCMC experiments.
### It contains the function for calculation of gradient and density of two gaussian mixture posterior. 
### And two additional functions to support the calculation. 


###============= twogaussian_par: supportive function 1==========
## calculation on the parameters of the multivariate normal, like the matrix inverse, so avoid repeating such calculation in sampling. 

## input: weight, sigma1, sigma2. The overall posterior distribution = weight*N(mu1,sigma1) + (1 - weight)*N(mu2, sigma2). 

## output: 
#   invsigma1, invsigma2: inverse of the covariance matrix. 
#   scale1, scale2: the scale in front of the exp of normal. 
#=======================
twogaussian_par <- function(weight, sigma1, sigma2)
{
  dimension = dim(sigma1)[1]
  detsigma1 = det(sigma1)   
  # determinant of the variance matrix 1
  detsigma2 = det(sigma2)   
  # determinant of the variance matrix 2
  scale1 = weight/sqrt((2*pi)^dimension*detsigma1) 
  # so that scale1 * exp(-matrix norm of (x - mu1)/2) = weighted density of normal1. 
  scale2 = (1 - weight)/sqrt((2*pi)^dimension*detsigma2)
  # so that scale2 * exp(-matrix norm of (x - mu2)/2) = weighted density of normal2. 
  
  invsigma1 = solve(sigma1) 
  # inverse of the variance matrix 1
  invsigma2 = solve(sigma2)
  # inverse of the variance matrix 1
  
  return(list("scale1" = scale1, "scale2" = scale2, "invsigma1" = invsigma1, "invsigma2" = invsigma2))
}




###============= calculate_density_grad: supportive function 2==========

## calculate_density_grad: a function calculate the posterior density and gradient of log density wrt a given x
## input:
#    ParResult:     The outputs of twogaussian_par. 
#    mu1, mu2: The overall posterior distribution = weight*N(mu1,sigma1) + (1 - weight)*N(mu2, sigma2). 
#    x: The current position of a sample.  

## output:
#    pd: The posterior density of x. 
#    gradU: The gradient of -log posterior. 

#======================
calculate_density_grad <- function(ParResult, mu1, mu2, x)
{
  scale1 = ParResult$scale1
  scale2 = ParResult$scale2
  invsigma1 = ParResult$invsigma1
  invsigma2 = ParResult$invsigma2
  
  xtomu1 = x - mu1
  xtomu2 = x - mu2
  
  weighted_pd1 = scale1*exp(-t(xtomu1) %*% invsigma1 %*% xtomu1/2) 
  # weight*N(mu1, sigma1)
  weighted_pd2 = scale2*exp(-t(xtomu2) %*% invsigma2 %*% xtomu2/2)
  # (1-weight)*N(mu2, sigma2)
  weighted_pd1 = as.numeric(weighted_pd1)
  weighted_pd2 = as.numeric(weighted_pd2)
  
  pd = weighted_pd1 + weighted_pd2
  # probability density at x. 
  gradU = ( weighted_pd1 * invsigma1 %*% xtomu1 + weighted_pd2 * invsigma2 %*% xtomu2)/pd
  # gradient of U at x. U = -log(probability density)
  
  return(list("pd" = pd, "gradU" = gradU))
}


### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
## parameter setting. 
dimension = 5
mu1 = rep(0,dimension)
mu2 = rep(5,dimension)
sigma1 = diag(1,dimension)
sigma2 = diag(1,dimension)
weight = 1/3
ParResult = twogaussian_par(weight, sigma1, sigma2)

# ===================================
## calculate the probability density and gradient for a position x.
## input: 
#    x:   positions. matrix of p*n. p = dimension of space, n = number of chains. 
## output: 
#    a data frame of two columns
#         $grad = gradient of the potential at x. 
#         $pd   = (posterior) probability density at x. 

## side note: potential = -log(probability density)
twogaussian <- function(x) 
{
  temp <- function(x)
  {
    temp = calculate_density_grad(ParResult, mu1, mu2, x)
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
# mu1 = c(1,1,0,0)
# mu2 = c(0,0,-1,-1)
# sigma1 = diag(1,4)
# sigma2 = diag(1,4)
# ParResult = twogaussian_par(1/2, sigma1, sigma2)
# ParResult$scale2
# ParResult$invsigma2
# 
# x = matrix(c(0,0,0,0,0,5,5,5,5,5),nrow = 5)
# result2 = twogaussian(result, mu1, mu2, c(1,1,0,0))
# result2$gradU



