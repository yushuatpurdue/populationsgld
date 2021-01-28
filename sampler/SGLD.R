## This is the one step SGLD sampler. 

## input: 
#      x:             current positions of the chains. 
# distribution_sgld:  target distribution. 
#epsilon:             learning rate. 
#    tau:             temperature

## output: 
#      positions after one step update of SGLD sampler. 

###=========== redundant block: should have realized in another file where this file is called ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/2gaussian.R")                   # source some potential
# distribution_sgld <-  twogaussian                   # pass it to a general name that is recognizable for function SGLD. 
# rm(twogaussian)                                     # remove the redundant source function. 

### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================

SGLD <- function(r, distribution_sgld, epsilon = 0.1,tau = 1){
  grad = distribution_sgld(r)$grad
  wn_length = length(r)
  wn = rnorm(wn_length)
  wn = array(wn, dim = dim(r))
  
  r_update = r - grad*epsilon +sqrt(2*epsilon*tau)*wn
  return(r_update)
}



##====================================
## =========== test ==================
##====================================

### test for SGLD applied on 
# r = array(c(1,-0.5,0.3,0.2,0.4),dim = c(1,5))
# epsilon = 0.1
# SGLD(r,epsilon)



