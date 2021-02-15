### This file is created by Yushu Huang for population SGMCMC experiments.
##  This is the one step snooker sampler. Main function named snooker, and a bunch of supportive functions. ##  Follows the Replica Exchange algorithm (Wei Deng 2020). 

##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: global parameter setups ================================
Ts                  = c(1,10)                    # temperature of the two chains. 
loglike_variance    = 0.7^2                        # tuning parameter: the variance of the approximated loglikelihood. In the mixture gaussian setting it cannot be estimated, so just fix and tune it. No need of F (correction factor) here consequently. 

### ========== call function: calculating potential & gradient ===============

source("potential/mixtureGaussian.R")  
potential = mixgaussian

# ==================================
## load the functions to be used when sampling. 

source("sampler/SGLD.R")                          # source 1 step SLGD as updating method. 

##===================================================
## =========== supportive function ==================
##===================================================

### ============ calculating the swapping rate ===========
# ========================================================
swaprate <- function(loglikelihood){
  # input:   loglikelihhod            = the loglikelihood of the current position of the two chains. length 2 vector. 
  difftemperature_term  = -diff(1/Ts)          # 1/t1 - 1/t2
  diffloglikelihood     = -diff(loglikelihood)
  swaprate              = exp(difftemperature_term * (diffloglikelihood - difftemperature_term*loglike_variance))
  return(swaprate)
}




