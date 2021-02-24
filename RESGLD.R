### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the function that runs replica exchange SGLD. 

## input: 
#     niteration         : number of iterations. scaler. 
#     learning rate      : learning rate for the two chains respectively. vector. 
#     Ts                 : temperatures for the two chains respectively. vector. 

## output: 
#     a population of paths. array of dimension = c(dim(state space), 2, niteration)

RESGLD <-function(Ts, niteration, learning_rate){
  
###=========== set path ================
  # setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  
###============= initialization: call functions and initial values=======================
  
  source("potential/mixtureGaussian.R")
  distribution_sgld <- mixgaussian               # target posterior distribution.  
  
  ## sample x. need to change this according to the target distribution!!
  # temp <- mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))      # initial states
  # x <<- t(temp)                                  # DEFINE AS GLOBAL! or snooker(x, Index) wouldn't be able to read it.
  temp <- runif(2*dimension, min = 0, max = 1)
  x <<- matrix(temp, nrow = dimension)
  print(x)
  
  source("sampler/replica.R")                    # source one-step replica exchange sampler. Note: SGLD sampler has got sourced in "sampler/replica.R", no need to source again. 
  
###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,2,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 2:niteration) {
    
      x = replica(x, learning_rate, Ts)
      output[,,t] = x
    
    # print(x)
    # print(t)
  }
  
  return(list("output" = output, "mu" = mu))
  
}



###=========================================
###=============== run function ============
###=========================================
Ts                  = c(1,10)                    # temperature of the two chains. 
learning_rate       = c(0.01,0.03)


output_RESGLD = RESGLD(Ts, niteration = 50000, learning_rate)

hist(as.numeric(output_RESGLD$output[1,2,]),breaks = seq(-7,12,0.1),main = "RESGLD")

plot(output_RESGLD$output[1,,], output_RESGLD$output[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "RESGLD 50000 iteration", xlab = "x", ylab = "y")
points(output_RESGLD$mu[1,], output_RESGLD$mu[2,], pch = 3,col = "green")

plot(output_RESGLD$output[1,1,], output_RESGLD$output[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "RESGLD: chain 1, 50000 iteration", xlab = "x", ylab = "y")
points(output_RESGLD$mu[1,], output_RESGLD$mu[2,], pch = 3,col = "green")

plot(output_RESGLD$output[1,2,], output_RESGLD$output[2,2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "RESGLD: chain 2, 50000 iteration", xlab = "x", ylab = "y")
points(output_RESGLD$mu[1,], output_RESGLD$mu[2,], pch = 3,col = "green")

write.csv(output_RESGLD$output, "./output_RESGLD.csv")
write.csv(output_RESGLD$mu, "./mu.csv")





