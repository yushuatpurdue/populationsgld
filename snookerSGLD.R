### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the function that runs population SGLD with snooker. 

## input: 
#     niteration   : number of iterations. 
#     nchain       : number of the chains in the population to run. 

## output: 
#     a population of paths. array of dimension = c(dim(state space), nchain, niteration)

###==============================================
###=============== !!!MAIN FUNCTION!!! ==========
###==============================================
library(MASS)

run_snookerSGLD <- function(Ts, niteration, nchain){
###=========== set path ================
  # setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  
###============= initialization: call functions and initial values=======================
  source("potential/mixtureGaussian.R")
  distribution_sgld <- mixgaussian               # target posterior distribution. 
  
  ## sample x. need to change this according to the target distribution!!
  # temp <- mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))      # initial states
  # x <<- t(temp)                                  # DEFINE AS GLOBAL! or snooker(x, Index) wouldn't be able to read it.
  temp <- runif(nchain*dimension, min = 0, max = 10)
  x <<- matrix(temp, nrow = dimension)
  
  source("sampler/snooker.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 1:niteration) {
    full_index = 1:nchain
    Index = sample.int(nchain, ceiling(nchain/2))        # sample 5 integers between 1 and 10 without replacement. 
    sgld_index = full_index[-Index]
    
    
    x[,Index] = snooker(x, Ts, Index)
    
    
    x[,sgld_index] = SGLD(x[,-Index], distribution_sgld, min(0.1,sigma^2))

    output[,,t] = x
    print(t)
  }
  
  return(list("output" = output, "mu" = mu))
  
}


###=========================================
###=============== run function ============
###=========================================

output_snooker = run_snookerSGLD(Ts, 5000,20)
# hist(as.numeric(output_snooker$output[2,,1:10]),breaks = seq(-7,12,0.5),main = "snooker SGLD")
write.csv(output_snooker$output, "./output_snooker.csv")
write.csv(output_snooker$mu, "./mu.csv")

###=========================================
###=============== unwrapped code ==========
###=========================================


# ###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# 
# nchain = 10
# niteration = 5000
# 
# ###============= initialization: call functions and initial values=======================
# source("potential/2gaussian.R")
# distribution_sgld <- twogaussian               # pass potential(r) to a general name that is recognizable for function SGLD.
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# source("sampler/snooker.R")
# 
# ###=============== to be wrapped into a function =====================
# 
# 
# output = array(NA,dim = c(2,10,niteration))
# output[,,1] = x
# 
# for (t in 1:niteration) {
#   full_index = 1:10
#   Index = sample.int(10,5)        # sample 5 integers between 1 and 10 without replacement.
#   sgld_index = full_index[-Index]
# 
#   # distribution_sgld <- snooker_scaler               # pass potential(r) to a general name that is recognizable for function SGLD.
# 
#   x[,Index] = snooker(x,Index)
# 
# 
#   x[,sgld_index] = SGLD(x[,-Index], distribution_sgld)
# 
#   output[,,t] = x
# }



