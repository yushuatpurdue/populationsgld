library(MASS)

Data = read.csv("dataset/datasetNLR1.csv")
source("potential/nonlinearReg.R")
# potential = nonlinearReg
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)
## note: need to specify batchs and nData. 

snookerSGLDbatch <- function(Ts, temperature_select, nepoch, nchain, learning_rate){
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
  # temp <- runif(nchain*dimension, min = 0, max = 1)
  # x <<- matrix(temp, ncol = dimension)
  theta <<- t(mvrnorm(nchain, mu = c(0,0), Sigma = diag(1,2)))
  # theta <<- t(matrix(c(20,10,20,10),ncol = 2))
  
  source("sampler/snookerBatch.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
  ###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,nepoch*nBatch))   # result array
  # output[,,1] = x                               # record initialization
  
  for (t in 1:nepoch) {
    for (k in 1:nBatch) {
      Batch = Data[which(Data$augment == k),]
      potential <<- function(theta){
        potential = stochPotential(theta, Batch, nData)
      }
      full_index = 1:nchain
      if (mutation_rate>0 && mutation_rate<1){
        Index = sample.int(nchain, ceiling(nchain*(1 - mutation_rate)))        # sample (1 - mutation_rate) portion of the chains to do crossover. 
        # print(Index)
        sgld_index = full_index[-Index]
        
        
        theta[,Index] = snooker_batch(theta, temperature_select, Index, tau = Ts)
        
        
        theta[,sgld_index] = SGLD(theta[,-Index], distribution_sgld = potential, epsilon = learning_rate[sgld_index], tau = Ts[sgld_index])
      }
      else if (mutation_rate == 1){
        theta = SGLD(theta, distribution_sgld = potential, epsilon = learning_rate, tau = Ts)
      }
      
      
      output[,,((t-1)*nBatch+k)] <- theta
      # print(theta)
    }
    
    
    
    # print(t)
  }
  
  return(list("output" = output))
  
}

## test
# Batch = Data[which(Data$augment == 1),]
# nData = 10000
# potential <- function(theta){
#   temp = nonlinearReg(theta, Batch, nData)
# }
# theta = mvrnorm(20, mu = c(0,0), Sigma = diag(1,2))
# nchain  = 20
# nepoch = 10
# nchain  = 10
Ts      = seq(1,5,4/(nchain-1))
# Ts      = rep(1,nchain)
mutation_rate = 0.6        # number of chains to go through SGLD instead of snooker.

learning_rate = seq(0.005,0.025,0.02/(nchain - 1))
output_snookerBatch = snookerSGLDbatch(Ts, temperature_select = 0.5, nepoch = 40, nchain = 10, learning_rate)
