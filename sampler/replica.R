### This file is created by Yushu Huang for population SGMCMC experiments.
##  This is the one step snooker sampler. Main function named snooker, and a bunch of supportive functions. ##  Follows the Replica Exchange algorithm (Wei Deng 2020). 

##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: global parameter setups ================================
loglike_variance    = 0.3^2                        # tuning parameter: the variance of the approximated loglikelihood. In the mixture gaussian setting it cannot be estimated, so just fix and tune it. No need of F (correction factor) here consequently. 


### ========== call function: calculating potential & gradient ===============

# source("potential/mixtureGaussian.R")  
# potential = mixgaussian

# ========================================================
## load the functions to be used when sampling. 

source("sampler/SGLD.R")                          # source 1 step SLGD as updating method. 

##========================================================
## ============= supportive function =====================
##========================================================

### ============ calculating the swapping rate ===========
# ========================================================
swaprate <- function(x, Ts){
  # input:   
  # x            = the current position of the two chains. dimension*2 matrix. 
  # Ts           = temperature of the two chains. length 2 vector. 
  # output:  swap rate of x. 
  
  # first calculate the gradient and pobability density of x. 
  pd_grad               = potential(x)
  likelihood            = pd_grad$pd
  loglikelihood         = log(likelihood)
  difftemperature_term  = -diff(1/Ts)          # 1/t1 - 1/t2
  diffloglikelihood     = -diff(t(loglikelihood))
  swaprate              = exp(difftemperature_term * (diffloglikelihood - difftemperature_term*loglike_variance))
  return(swaprate)
}


##===================================================
## =========== !!!MAIN FUNCTION!!! ==================
##===================================================
## This is the one step replica exchange sampler. 

## input: 
#     x:                       current positions of the two chains dim(x) = [dimension, 2]
#     learning rate:           the learning rate for each chain. (a scaler or) a length 2 vector. 
#     Ts:                      temperature of the two chains. (a scaler or) a length 2 vector. 

## output: 
#     those states of x that got updated.
#====================

replica <-function(x, learning_rate = 0.1, Ts){
  
  # 1) one step SGLD
  x      = SGLD(x, potential, learning_rate, Ts)
  
  # 2) Decided to swap or not. 
  swap_x = swaprate(x, Ts)          # swapping rate. 
  u      = runif(1)              # generate an U[0,1] random number
  if (u < swap_x){              # swap if u < swapping rate. 
    temp = x[,1]
    x[,1] <- x[,2]
    x[,2] <-temp
  }
  
  return(x)
}





##====================================
## =========== test ==================
##====================================
# dimension = 2
# x = mvrnorm(2, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# learning_rate = c(0.01,0.03)
# Ts = c(1,10)                    # temperature of the two chains.
# replica(x,learning_rate, Ts)
