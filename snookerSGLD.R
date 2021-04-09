### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the function that runs population SGLD with snooker. 

## input: 
#     Ts                 : temperature ladder. vector. 
#     temperature_select : select temperature in snooker. scaler. 
#     niteration         : number of iterations. scaler. 
#     nchain             : number of the chains in the population to run. scaler. 

## output: 
#     a population of paths. array of dimension = c(dim(state space), nchain, niteration)

###==============================================
###=============== !!!MAIN FUNCTION!!! ==========
###==============================================
library(MASS)

snookerSGLD <- function(Ts, temperature_select, niteration, nchain, learning_rate){
###=========== set path ================
  # setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  
###============= initialization: call functions and initial values=======================
  source("potential/mixtureGaussian.R")
  distribution_sgld <- mixgaussian               # target posterior distribution.  
  ## sample x. need to change this according to the target distribution!!
  # temp <- mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))      # initial states
  # x <<- t(temp)                                  # DEFINE AS GLOBAL! or snooker(x, Index) wouldn't be able to read it.
  temp <- runif(nchain*dimension, min = 0, max = 1)
  x <<- matrix(temp, nrow = dimension)
  
  source("sampler/snooker.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 2:niteration) {
    full_index = 1:nchain
    if (mutation_rate>0 && mutation_rate<1){
      Index = sample.int(nchain, ceiling(nchain*(1 - mutation_rate)))        # sample (1 - mutation_rate) portion of the chains to do crossover. 
      # print(Index)
      sgld_index = full_index[-Index]
      
      
      x[,Index] = snooker(x, temperature_select, Index)
      
      
      x[,sgld_index] = SGLD(x[,-Index], distribution_sgld, epsilon = learning_rate, tau = Ts[sgld_index])
      
      output[,,t] = x
    }
    else if (mutation_rate == 1){
      x = SGLD(x, distribution_sgld, epsilon = learning_rate, tau = Ts)
      output[,,t] = x
    }
   
    # print(x)
    # print(t)
  }
  
  return(list("output" = output, "mu" = mu))
  
}


###=========================================
###=============== run function ============
###=========================================
nchain  = 20
Ts      = seq(1,5,4/(nchain-1))
mutation_rate = 0.8        # number of chains to go through SGLD instead of snooker.
learning_rate = seq(0.01,0.03,0.02/(nchain*mutation_rate - 1))


##================ run option1: evaluate the mean & var estimation ==================
# nrun    = 4                  # run nrun independent times. 
# mu_est_EMC = array(NA, dim = c(2, nrun))
# var_est_EMC = array(NA, dim = c(2,2,nrun))
# for (i in 1:nrun) {
#   output_snooker = snookerSGLD(Ts, temperature_select = 0.5, 200000, nchain, learning_rate)
#   # hist(as.numeric(output_snooker$output[2,,1:10]),breaks = seq(-7,12,0.5),main = "snooker SGLD")
#   record = output_snooker$output[,1,]
#   
#   mu_est_EMC[,i] = rowMeans(record)
#   var_est_EMC[,,i] = cov(t(record))
# }
# 
# print(mu_est_EMC)
# print(var_est_EMC)
# 
# jobid = Sys.getenv("SLURM_JOBID")
# write.csv(mu_est_EMC, paste0("./mu_est_EMC",jobid,".csv"))
# write.csv(var_est_EMC, paste0("./var_est_EMC",jobid,".csv"))


#
##==================== run option2: run chain and plot =================
 output_snooker = snookerSGLD(Ts, temperature_select = 0.5, 20000, nchain, learning_rate)
#
## plot(output_snooker$output[1,,], output_snooker$output[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker 20 chains 200000 iteration", xlab = "x", ylab = "y")
## points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
#
plot(output_snooker$output[1,1,], output_snooker$output[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 20000 iteration", xlab = "x", ylab = "y")
 points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
## 
## plot(output_snooker$output[1,1,1:50000], output_snooker$output[2,1,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 50000 iteration", xlab = "x", ylab = "y")
## points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
#
write.csv(output_snooker$output[,1,], "./output_snooker.csv")



##===========================
## run just SGLD:
# nchain = 1
# Ts = 0.8
# mutation_rate = 1
# learning_rate = 0.02
# output_SGLD = snookerSGLD(Ts, temperature_select = 0.5, 200000, nchain, learning_rate)
# 
# output = array(output_SGLD$output,dim = c(2,200000))
# output = as.matrix(output)
# 
# plot(output[1,1:50000], output[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGLD 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output[1,], output[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGLD 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# write.csv(output, "./variousSGLD/output_SGLD.csv")

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




