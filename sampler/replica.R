### This file is created by Yushu Huang for population SGMCMC experiments.
##  This is the one step snooker sampler. Main function named snooker, and a bunch of supportive functions. ##  Follows the Replica Exchange algorithm (Wei Deng 2020). 

##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: parameter setups ===================================================
Ts      = seq(1,5,4/(nchain-1))            # selection temperature for sampling xj

### ========== call function: calculating potential & gradient ===============

source("potential/mixtureGaussian.R")  
epsilon_r = min(0.1, sigma^2)                        # learning rate for sgld in sampling r. 
potential = mixgaussian

# ==================================
## load the functions to be used when sampling. 

source("potential/snooker_scaler.R")              # source the potential of r. 

source("sampler/SGLD.R")                          # source 1 step SLGD as updating method. 


