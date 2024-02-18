### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation of gradient and density of Gaussian mixture posterior of arbitrary number of Gaussians. 
### And one additional function to support the calculation. 

##### Note that the mixture Gaussians have equivariance of a diagonal form. So this file does not apply to the most general case. 



###============= calculate_density_grad: supportive function ==========

## calculate_density_grad: a function calculate the posterior density and gradient of - log density wrt a given x
## input:
#    mu:      Matrix of size (d, m). Column i = mean vector of i-th Gaussian. m = number of Gaussian in the model. 
#    sigma:   Scaler. Each Gaussian has the same variance of diag(sigma^2,d)
#    x:       The current position of a sample.  
#    weight:  (p1, p2, p3,..., pm), a m*1 vector of weights sums up to 1, each entry a weight. 

## output:
#    $pd: The posterior density at x. 
#    $gradU: The gradient of -log posterior. 

#======================
calculate_density_grad <- function(mu, sigma, x, weight)
{
  dimension = dim(mu)[1]
  n_mixture = dim(mu)[2]
  mixture_scale = weight/(sqrt(2*pi)*sigma^dimension) 
  # so that scale * exp(-matrix norm of (x - mu)/2) = a vector of weighted density of each Gaussian.  
  
  xtomu = x - mu  
  # matrix of dimension * #Gaussian
  variance = sigma^2
  l2dist = colSums(xtomu^2)
  # a vector of length m = #Gaussians. 
  
  # left_xtomu = array(xtomu, dim = c(length(xtomu),1))
  # Sigma_wrap = bdiag(apply(invSigma, n_mixture, bdiag))
  # right_xtomu = bdiag(split(xtomu, rep(1:n_mixture, each = dimension)))
  
  Gaussian_dist = exp(-l2dist/variance/2)
  # the values of the Gaussian kernel within each of the Gaussians. 
  
  weighted_pd = mixture_scale * Gaussian_dist 
  weighted_pd = as.numeric(weighted_pd)
  # weighted probability densities of individual Gaussians in the mixture model. 
  pd = sum(weighted_pd)
  # probability density at x. 

  gradU = xtomu%*%weighted_pd/variance/pd
  # gradient of U at x. U = -log(probability density)
  
  return(list("pd" = pd, "gradU" = gradU))
}


### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
## parameter setting. 
epsilon_r = 0.003  

dimension = 2
n_mixture = 20

weight = rep(1/n_mixture, n_mixture)
sigma = 0.1
mu = runif(dimension*n_mixture,0,10)
mu = matrix(mu, nrow = dimension)
# mu = t(matrix(c(2.418759, 6.618219, 6.691523, 6.904582, 2.196740, 0.5989215, 0.1847932, 5.264697, 4.879505, 8.004288, 6.077461, 3.801844, 5.177598, 2.469589, 9.876724, 2.1504432, 3.2653869, 2.133480, 7.884069, 4.334636),ncol = dimension))

# mu = t(matrix(c(2.547309, 1.005194, 6.898408, 9.168438, 7.154848, 9.028120, 2.140865, 6.635207, 9.1268483, 0.03306429, 9.616817, 1.475194, 1.15205, 8.724225, 3.652639, 8.746123, 7.246912, 4.365018, 6.467183, 2.876139,
#                 2.547301, 1.910528, 8.405147, 5.954380, 3.854593, 9.139374, 8.745735, 7.449630, 0.9170815, 5.54110136, 6.089943, 9.411089, 1.08722, 4.368250, 8.852135, 2.893911, 2.387444, 6.759031, 6.694622, 6.853250), ncol = dimension))    # ESGLD experiment

mu = t(matrix(c(0.09530785, 9.98696125, 9.07519804, 9.51948259, 2.88687893, 1.96513550, 2.74847636, 9.57410327, 6.14152122, 3.49452700, 9.44337555, 4.22920353, 2.74555399, 8.85162807, 3.81950937, 4.46274869, 8.41448503, 5.33006326, 8.09606520, 6.94301054,
                1.1296385, 8.2583106, 8.9988781, 3.7418591, 5.0763131, 2.0828439, 0.5741568, 0.2035272, 9.3144791, 3.4301498, 6.9505852, 0.2078327, 6.8422350, 1.3386339, 9.0890026, 9.2310132, 4.9809740, 1.6721373, 4.0875138, 1.9924000), ncol = dimension))   # SPOS experiment
mu = matrix(c(1,5,5,1,7,7), nrow = dimension, ncol = n_mixture)


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
    temp = calculate_density_grad(mu, sigma, x, weight)
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
# weight = c(1/2, 1/2, 1/3)
# Sigma = array(NA, dim = c(2,2,3))
# mu = matrix(NA, nrow = 2, ncol = 3)
# for (i in 1:3) {
#   Sigma[,,i] = diag(i,2)
#   mu[,i] = c(i,i)
# }

# x = t(mvrnorm(1, mu = rep(0,2), Sigma = diag(1,2)))
# 
# ParResult = mixgaussian_par(weight, Sigma)
# ParResult$mixture_scale
# ParResult$invSigma
# 
# # x = matrix(c(1,1,0,0, 0,0,1,1),nrow = 4)
# result = mixgaussian(x)
# result$grad



