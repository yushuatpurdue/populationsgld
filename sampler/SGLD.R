### This file is created by Yushu Huang for EMC-SGLD experiments.
### This is the one step SGLD sampler. 
### An additional supportive function is defined to deal with the gradient blowup. 


###=========== redundant block: should have realized in another file where this file is called ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/2gaussian.R")                   # source some potential
# distribution_sgld <-  twogaussian                   # pass it to a general name that is recognizable for function SGLD. 
# rm(twogaussian)                                     # remove the redundant source function. 

### ====================================================
### ======== supportive function: rowsum_mu=============
### ====================================================
## Idea: 
# 1. When a gradient vanishes in a not-too-sparse space, it must be too far away from the centers region, meaning that we need to find the direction of the centers region. 
# 2. When a gradient vanishes in a too-sparse space, it could either be too far away from where the centers are, or it could be amid the centers while too far from any of them, in this case we need to keep moving to any direction until the gradient escapes from being NA. 
# 3. Once moving into the correct direction the sample would be kicked into the centers region, and should be playing around quite well that the vanishing gradient doesn't happen again. 

## Solution: if gradient = NA, then calculate x - mu and sum it up by directions, i.e. rowSums(x - mu). If it's positive in a direction then x needs to be dragged towards negative, so set gradU = 200, else if it's negative then set gradU = -200. This is a bit similar to the reLU idea. 

### The function "rowsum_mu" is to calculate x - mu and sum it up by directions. 
## parameter needed: 
#     mu:          the center of the gaussians in the mixture gaussian model. 
## input: 
#     x:           a position, length d vector. 
## output: 
#     length d vector, each entry the sum of difference of the position and mu in that dimension. 

if (exists("mu")){
  rowsum_mu <- function(x){
    return(rowSums(x - mu))
  }
}




### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================

## input: 
#      r:             current positions of the chains. A d*nchains matrix.  
# distribution_sgld:  target distribution function, returns density and SGLD gradient. In "./potential/". 
# epsilon:            learning rate. 
#    tau:             temperature

## output: 
#      positions after one step update of SGLD sampler. 

SGLD <- function(r, distribution_sgld, epsilon = 0.1,tau = 1){
  grad = distribution_sgld(r)$grad      # compute the gradient. 
  
  rowsum_mu <- function(x){
    return(rowSums(x - mu))
  }
  
  if (sum(is.na(grad)) > 0){            # if there exists any entry being NA in grad:
    if (exists("mu")){
      NAdirection = sign(apply(r, 2, rowsum_mu))                   #  returns the direction for correction
      blowup_entry = NAdirection*(is.na(grad))         # (1-is.na(grad)) returns which entries blowed up. sign(rowSums(x - mu)) 
      
      grad[is.na(grad)] <- 0
      grad = 200*blowup_entry + grad
    }
    else {
      stop("gradient blows up")
    }
  }
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



