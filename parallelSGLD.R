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

library(MASS)

run_parallelSGLD <- function(niteration, nchain, epsilon = 0.1, tau = 1){
  ###=========== set path ================
  # setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  
  ###============= initialization: call functions and initial values=======================
  source("potential/2gaussian.R")
  source("sampler/SGLD.R")
  distribution_sgld <- twogaussian               # target posterior distribution. 
  temp <- matrix(mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension)), nrow = nchain)      # initial states
  x <<- t(temp)                                  # DEFINE AS GLOBAL! or SGLD(x, ...) wouldn't be able to read it. 
  
  
  ###=============== run parallel SGLD =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  output[,,1] = x                               # record initialization
  
  for (t in 1:niteration) {
    
    x = SGLD(x, distribution_sgld,epsilon = 0.1, tau = 1)
    
    output[,,t] = x
  }
  
  return(output)
  
}

###=========================================
###=============== run function ============
###=========================================




nrun    = 3                  # run nrun independent times.

for (i in 1:nrun) {
  start_time <- Sys.time()
  output_parallel = run_parallelSGLD(1.3e7,1,epsilon = 0.05)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))
  # write.csv(output_parallel, paste0("./output_parallel_",i, ".csv"))
  pdf(file=paste0("2gaussianSGLDplot_",i,".pdf"))
  for (j in 1:10) {
    bins = hist(as.numeric(output_parallel[1,1,1:(1.3e6*j)]), breaks = seq(-8,13,0.5), freq=FALSE, main = paste0("SGLD ", 1.3e6*j, " iterations"), xlab = "axis 1")
    
  }
  thining = seq(1.3e4,1.3e7,1.3e4)
  time = seq(19/1000,19,19/1000)
  plot(time, output_parallel[1,1,thining], ylim = c(-2,8), type = "l", xlab = "CPU time (mins)", ylab = "position", main = "SGLD sample path")
  
  dev.off() 
}




