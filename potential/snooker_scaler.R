### This file is created by Yushu Huang for population SGMCMC experiments.
### It contains the function for calculation of gradient and density of the distribution that snooker scaler r should be sampled from. 
### global variable used in this file: potential. Defined in "sampler/snooker.R". 

## input: 
#    r:           snooker scaler. array of dim = c(1, #{x_j})

## output: 
#    a data frame of two columns
#         $grad = gradient of the potential at x. 
#         $pd   = (posterior) probability density at x. 

###=========== redundant block: should have realized in another file where this file is called ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/2gaussian.R") 


### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
snooker_scaler <- function(r, distr_grad, paired_x, direction){
  
  ## 
  
  dimension = dim(paired_x)[1]
  nsample = dim(paired_x)[2]
  
  r_matrix = t(array(r, dim = c(length(r), dimension)))
  position = paired_x + r_matrix*direction                     # xj + r*e
  
  temp_potential = distr_grad(position)                # compute pi(xj + r*e) and \nabla U(xj + r*e)
  pd = temp_potential$pd                               # No temperature ladder here
  grad = temp_potential$grad                           # No temperature ladder here
  # hess = temp_potential$hess 
  
  r_to_dm1 = abs(r)^(dimension - 1)    
  pd_r = r_to_dm1 *pd                                   # compute f(r) UP TO A CONSTANT!!
  numerator = colSums(direction*grad)
  grad_r = -(dimension - 1)/r + numerator                       # compute \nabla_r (-log(f(r) )) 
  # 
  # term1_hess = (dimension - 1)/r^2
  # hess_r = rep(NA,nsample)
  # for (i in 1:nsample) {
  #   h = term1_hess[i] + t(direction[,i])%*%hess[,,i]%*%direction[,i]
  #   hess_r[i] = h
  # }
  
  result = list("grad" = grad_r,"pd" = pd_r)
  return(result)
}


##====================================
## =========== test ==================
##====================================

# r = rnorm(length(Index))
# r = array(r,dim = c(1,length(r)))
# nchain = 10
# dimension = 5
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# 
# paired_x = x[,paired_Index]
# direction = direction_e
# snooker_scaler(r, paired_x = paired_x, direction = direction)
# rm(paired_x)
# rm(direction)

# distribution plot
r = array(seq(-10,10,0.05),dim = c(1,401))
r_matrix = t(array(r, dim = c(length(r), dimension)))
position = paired_x[,3] + r_matrix*direction [,3]
temp_potential = twogaussian(position)                # compute pi(xj + r*e) and \nabla U(xj + r*e)
pd = temp_potential$pd
grad = temp_potential$grad
r_to_dm1 = abs(r)^(dimension - 1)
pd_r = r_to_dm1 *pd                                   # compute f(r) UP TO A CONSTANT!!
numerator = colSums(direction[1]*grad)
grad_r = -(dimension - 1)/r + numerator

plot(r, pd_r, type = "l")
plot(r[2:401],-(log(pd_r)[2:401] - log(pd_r)[1:400])/0.05, type = "l")
plot(r, grad_r,type = "l")
plot(r,-(dimension - 1)*sign(r)/r,type = "l")
# 
# rm(r_matrix)
# rm(direction)
# rm(paired_x)
# rm(position)
