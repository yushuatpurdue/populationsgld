### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the function that runs basic SGLD parallelly. Serves as baseline 

## input: 
#     niteration   : number of iterations. 
#     nchain       : number of the chains in the population to run. 

## output: 
#     a population of paths. array of dimension = c(dim(state space), nchain, niteration)

###==============================================
###=============== !!!MAIN FUNCTION!!! ==========
###==============================================

run_parallelSGLD <- function(niteration, nchain){
  ###=========== set path ================
  setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  
  ###============= initialization: call functions and initial values=======================
  source("potential/2gaussian.R")
  distribution_sgld <- twogaussian               # target posterior distribution. 
  temp <- mvrnorm(nchain, mu = rep(2.5,dimension), Sigma = diag(1,dimension))      # initial states
  x <<- t(temp)                                  # DEFINE AS GLOBAL! or snooker(x, Index) wouldn't be able to read it. 
  
  source("sampler/snooker.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
  ###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 1:niteration) {
    
    x = SGLD(x, distribution_sgld)
    
    output[,,t] = x
  }
  
  return(output)
  
}

###=========================================
###=============== run function ============
###=========================================

output_parallel = run_parallelSGLD(20000,10)
hist(output_parallel[1,,],breaks = seq(-7,12,0.5), main = "parallel SGLD")


