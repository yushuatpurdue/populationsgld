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

source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian               # target posterior distribution.
potential = mixgaussian

snookerSGLD <- function(Ts, temperature_select, niteration, nchain, learning_rate){
###=========== set path ================
  # if(getwd() != "/Users/yushuhuang/Documents/pu/2020fall/populationsgld"){
  #   setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  # }
  
###============= initialization: call functions and initial values=======================
  ## source the potentials to be used for sampling, then define the distribution accordingly. 
  
  # source("potential/2gaussian.R")
  # distribution_sgld <- twogaussian               # target posterior distribution.  
  
  ## sample initial position x. need to change this according to the target distribution!!
  # temp <- mvrnorm(nchain, mu = rep(2.5,dimension), Sigma = diag(25,dimension))      # initial states
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
      
      
      x[,Index] = snooker(x, temperature_select, Index, tau = Ts)
      
      
      x[,sgld_index] = SGLD(x[,-Index], distribution_sgld, epsilon = learning_rate[sgld_index], tau = Ts[sgld_index])
      
      output[,,t] = x
    }
    else if (mutation_rate == 1){
      x = SGLD(x, distribution_sgld, epsilon = learning_rate, tau = Ts)
      output[,,t] = x
    }
   
    # print(x)
    # print(t)
  }
  
  return(list("output" = output))
  
}


###=========================================
###=============== run function ============
###=========================================
nchain  = 20
# nchain  = 10
Ts      = seq(1,5,4/(nchain-1))
# Ts      = rep(1,nchain)
mutation_rate = 0.90        # number of chains to go through SGLD instead of snooker.
print(paste0("mutation rate = ", mutation_rate))
# mutation_rate = 0.5        # number of chains to go through SGLD instead of snooker.
learning_rate = seq(0.01,0.03,0.02/(nchain - 1))
# learning_rate = rep(0.05,nchain)


##================ run option1: evaluate the mean & var estimation ==================
nrun    = 3                  # run nrun independent times.
mu_est_EMC = array(NA, dim = c(2, nrun))
var_est_EMC = array(NA, dim = c(2,2,nrun))

for (i in 1:nrun) {
  start_time <- Sys.time()
  output_snooker = snookerSGLD(Ts, temperature_select = 0.5, 20000, nchain, learning_rate)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))
  # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
  # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))
  
  
  record = output_snooker$output[,1,]

  mu_est_EMC[,i] = rowMeans(record)
  var_est_EMC[,,i] = cov(t(record))

  ## plot
#   pdf(file=paste0("2gaussianEMCplot_",i,".pdf"))
#   for (j in 1:20) {
#     hist(as.numeric(output_snooker$output[1,1,1:(100*j)]),breaks = seq(-8,13,0.5), freq=FALSE, main = paste0("EMC-SGLD ", 100*j, " iterations"), xlab = "axis 1")
#   }
#   plot(seq(22/2000,22,22/2000),output_snooker$output[1,1,0:2000], type = "l", ylim = c(-2,8), main = "EMC-SGLD sample path", xlab = "CPU time (mins)", ylab = "position")
#   
#   dev.off() 
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
# output_snooker = snookerSGLD(Ts, temperature_select = 0.5, 50000, nchain, learning_rate)
# # 
  pdf(file=paste0("mixGaussianEMCplot_",i,".pdf"))
    for (j in 1:20) {
      plot(output_snooker$output[1,1,1:(1000*j)], output_snooker$output[2,1,1:(1000*j)], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("snooker ",1000*j, " iteration"), xlab = "x", ylab = "y")
      points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
    
    }
  plot(output_snooker$output[1,1,0:20000],output_snooker$output[2,1,0:20000], type = "l", xlim = c(0,10), ylim = c(0,10), main = "EMC-SGLD sample path", xlab = "x", ylab = "y")
  
  dev.off() 
  print(paste("mu estimate = ", mu_est_EMC[,i]))
  print(paste("var estimate = ", var_est_EMC[,,i]))
}

# # plot(output_snooker$output[1,,], output_snooker$output[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker 20 chains 200000 iteration", xlab = "x", ylab = "y")
# # points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# # 
# plot(output_snooker$output[1,1,], output_snooker$output[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 50000 iteration", xlab = "x", ylab = "y")
#  points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# 
# plot(output_snooker$output[1,1,1:50000], output_snooker$output[2,1,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 50000 iteration", xlab = "x", ylab = "y")
# points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# 
# write.csv(output_snooker$output[,1,], "./output_snooker.csv")



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




