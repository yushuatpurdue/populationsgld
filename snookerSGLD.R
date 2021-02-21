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
  temp <- runif(nchain*dimension, min = 0, max = 10)
  x <<- matrix(temp, nrow = dimension)
  
  source("sampler/snooker.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 2:niteration) {
    full_index = 1:nchain
    if (mutation_rate>0 && mutation_rate<1){
      Index = sample.int(nchain, ceiling(nchain*(1 - mutation_rate)))        # sample (1 - mutation_rate) portion of the chains to do crossover. 
      print(Index)
      sgld_index = full_index[-Index]
      
      
      x[,Index] = snooker(x, temperature_select, Index)
      
      
      x[,sgld_index] = SGLD(x[,-Index], distribution_sgld, epsilon = learning_rate, tau = Ts[sgld_index])
      
      output[,,t] = x
    }
    else if (mutation_rate == 1){
      x = SGLD(x, distribution_sgld, epsilon = learning_rate, tau = Ts)
      output[,,t] = x
    }
   
    print(x)
    print(t)
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

output_snooker = snookerSGLD(Ts, temperature_select = 0.5, 10000, nchain, learning_rate)
# hist(as.numeric(output_snooker$output[2,,1:10]),breaks = seq(-7,12,0.5),main = "snooker SGLD")

plot(output_snooker$output[1,,], output_snooker$output[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker 20 chains 10000 iteration", xlab = "x", ylab = "y")
points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")

plot(output_snooker$output[1,1,], output_snooker$output[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 10000 iteration", xlab = "x", ylab = "y")
points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")

plot(output_snooker$output[1,5,], output_snooker$output[2,5,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 5, 10000 iteration", xlab = "x", ylab = "y")
points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")

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




# mu = matrix(c(0.526529485359788,6.22532256180421,0.104487277567387,2.02863411745057,7.69018175778911,3.42775387922302,4.20330020366237,3.67638011928648,2.2684468049556,6.03638118831441,0.0909082987345755,5.18603075295687,0.792325402144343,8.7226249021478,5.49441439798102,5.4901418578811,9.32448693783954,4.26261536078528,4.35691876569763,1.13028466468677),nrow = 2)


# x = matrix(c(3.499208, 3.601590, 3.560989, 3.282698, 3.378614, 3.436766, 3.523944, 3.460612, 3.487643, 3.488807, 3.442904, 7.225471, 3.562079, 3.549000, 3.435151, 3.464369, 3.650594, 3.398530, 3.452040, 3.360090, 9.738214, 9.914768, 9.851477, 9.791644, 9.704960, 9.806821, 9.826037, 9.690990, 9.787247, 9.763583, 9.894972, 9.549223, 9.617040, 9.882257, 9.698107, 9.892362, 9.827886, 9.807196, 9.609892, 9.875275), nrow = 2)
# Index = c(17,  7, 10, 14, 15, 19,  1,  8, 16,  3)
