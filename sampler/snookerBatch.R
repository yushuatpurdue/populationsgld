### This file is created by Yushu Huang for population SGMCMC experiments.
##  This is the one step snooker sampler. Main function named snookerBatch, and a bunch of supportive functions. 
##  This sampler is designed for the optimization of parameter on MINIBATCH dataset. 
##  Follows the EMC algorithm (Faming Liang 2001) and reconstruct it under the SGLD setting. 


##=======================================================
## =========== structural introduction ==================
##=======================================================
## THIS CODE CONSISITS OF 3 PARTS. 

## part 1: initialization. Many of this are inheritated from EMCSGLDBatch.R

## part 2: supportive functions. 
# Boltzmann        : calculating the boltzmann weights of the current chains, for the subsequent sampling. 
# sample_Boltzmann : sample the paired chains xj according to the Boltzmann weight. 
# l2               : calculate the l2 distance, to normalize xj - xi
# unit_direction   : calculate normalized direction of xj - xi
# sample_r         : sample the snooker scaler. 

## part 3: main function: snooker
# for a population of chains, uniformly draw a chain xi (this step is not included in this file), and pair it with xj chain drawn under Boltzmann distribution. 
# calculate the unit direction, 
# sample the snooker scaler with batches. 
# updated xi = xj + r*e 


##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: parameter setups ===================================================
nchain  = dim(theta)[2]                    # number of chains
# temperature_select = 0.5             # selection temperature for sampling xj
# Ts      = seq(1,5,4/(nchain-1))      # temperature ladder. Don't use this in snooker: inheritate from EMCSGLDBatch.R
### ========== call function: calculating potential & gradient ===============

# source("potential/mixtureGaussian.R")                # source the potential & gradient calculator.
epsilon_r = 0.003                        # learning rate for sgld in sampling r.
# stochPotential = nonlinearReg
epochR = 20
# Data = read.csv("dataset/datasetNLR1.csv")
# nData = dim(Data)[1]
# nBatch = max(Data$augment)

# source("potential/2gaussian.R")                # source the potential & gradient calculator.
# epsilon_r = 0.05  
# potential = twogaussian
## potential will be specified in snookerSGLD. 
# ==================================
## load the functions to be used when sampling. 

source("potential/snooker_scaler.R")              # source the potential of r. 

source("sampler/SGLD.R")                          # source 1 step SLGD as updating method. 

##====================================================
## =========== supportive functions ==================
##====================================================

### ============ Second: find the chains to pair with x ===========
# ==================================
## function calculate the boltzmann weight. 

Boltzmann <- function(x, temperature_select){
  # input:  x                      = vector of positions (again vector), so overall matrix
  #         temperature_select     = temperature selected for sampling, scaler
  # output: Boltzmann(x)           = their boltzmann weight = exp(-u(x)/temperature_select), vector
  pd_grad = potential(x)
  Boltzmann = (pd_grad$pd)^(1/temperature_select)
  return(Boltzmann)
}

# ==================================
## sample xj by boltzmann weight. 
sample_Boltzmann <- function(x, temperature_select, index){
  # input:  x                   = vector of positions (again vector), so overall matrix
  #         temperature_select  = temperature selected for sampling, scaler
  #         index               = the indices that should be take away before sampling
  # output: sampleid            = the indices of the sampled chain in {x} to snooker with xi
  sample_Boltzmann = c()
  for (i in index) 
  {
    weight = Boltzmann(x, temperature_select)     # calculate the boltzmann weights for everyone in {x}
    weight_net = weight[-i]   # take away xi
    weight_normalized = weight_net/sum(weight_net)   # normalize the weight of {x\xi} set
    samples = rmultinom(1,1,weight_normalized)   # sample according to the weight
    sampleid = which(samples>0)      # which one among the {x\xi} set is sampled?
    addon = as.numeric(sampleid>=i)
    sampleid = addon+sampleid        # what is the index of the sampled one in {x}? 
    sample_Boltzmann   = rbind(sample_Boltzmann, sampleid)
  }
  sample_Boltzmann = as.vector(sample_Boltzmann)
  
  # result = list(samples,sampleid,addon,sample_Boltzmann)
  return(sample_Boltzmann)
}



### ============ third: calculate the unit direction vector e ===========
# ==================================
## define a function to calculate the l2 norm of a vector
l2 <- function(vec){
  vec2 = vec^2
  l2   = sqrt(sum(vec2))
  return(l2)
}

# ==================================
## function to calculate e: the unit updating direction. 
unit_direction <- function(x, index, pair){
  diff_pair = x[,pair] - x[,index]
  e_norm = apply(diff_pair,2,l2)
  e_norm = rep(1,dimension)%*%t(e_norm)
  unit_direction = diff_pair/e_norm 
  return(unit_direction)
}


### ============ fourth: sample r ===========

# ==================================
## function: r sampler. 
sample_r <- function(epsilon, epochR = 20, temperature, T_annealing = 0.998, e_annealing = 0.999, Index,paired_Index, direction_e){
  # initializing r.  
  r = rnorm(length(Index))
  r = array(r,dim = c(1,length(r)))
  # initializing temperature and stepsize (for annealing), assuming temperature >= 1
  # temperature0 = temperature - 1 + 1/(T_annealing)^(0.5*epochR*nBatch)
  # anneal_temperature = temperature0
  # epsilon_r0 = epsilon_r/(e_annealing)^(0.5*epochR*nBatch)
  # anneal_epsilon_r = epsilon_r0
  # R = matrix(ncol = 5)
  # R[1,] = r
  
  # Annealing: half of total runs. 
  # for (t in 1:(0.5*epochR)){
  #   for (k in 1:nBatch){
  #     Batch = Data[which(Data$augment == k),]
  #     temp <- function(theta){
  #       temp = stochPotential(theta, Batch, nData)
  #     }
  #     distribution_snookerscaler <- function(r){
  #       distribution_snookerscaler = snooker_scaler(r,temp,theta[,paired_Index],direction_e)
  #     }
  #     anneal_epsilon_r = anneal_epsilon_r*e_annealing
  #     dist_to0 = min(abs(r),epsilon = anneal_epsilon_r)               # the learning rate needs to be small near 0 to avoid blowing up there.  
  #     anneal_temperature = (anneal_temperature+1-temperature)*T_annealing + temperature-1
  #     r = SGLD(r, distribution_snookerscaler, dist_to0, tau = anneal_temperature)
  #   }
  # }
  for (t in 1:(epochR)){
    for (k in 1:nBatch){
      Batch = Data[which(Data$augment == k),]
      temp <- function(theta){
        temp = nonlinearReg(theta, Batch, nData)
      }
      distribution_snookerscaler <- function(r){
        distribution_snookerscaler = snooker_scaler(r,temp,theta[,paired_Index],direction_e)
      }
      dist_to0 = min(abs(r),epsilon = epsilon_r)               # the learning rate needs to be small near 0 to avoid blowing up there.  
      r = SGLD(r, distribution_snookerscaler, dist_to0, tau = temperature)
    }
  }
  return(r)
}

# test for the convergence of sampling r.
# r0 = rnorm(length(Index))
# r0 = array(r0,dim = c(1,length(r0)))
# R = matrix(ncol = length(r0))
# R[1,] = r0
# r = r0
# epsilon_r = 0.01
# # T_annealing = 1
# epochR = 5
# # e_annealing = 1
# temperature = 1
# # temperature0 = temperature - 1 + 1/(T_annealing)^(0.5*epochR*nBatch)
# # anneal_temperature = temperature0
# # epsilon_r0 = epsilon_r/(e_annealing)^(0.5*epochR*nBatch)
# # anneal_epsilon_r = epsilon_r0
# # R = matrix(ncol = 5)
# # R[1,] = r
# 
# # Annealing: half of total runs.
# # for (t in 1:(0.5*epochR)){
# #   for (k in 1:nBatch){
# #     Batch = Data[which(Data$augment == k),]
# #     temp <- function(theta){
# #       temp = nonlinearReg(theta, Batch, nData)
# #     }
# #     distribution_snookerscaler <- function(r){
# #       distribution_snookerscaler = snooker_scaler(r,temp,theta[,paired_Index],direction_e)
# #     }
# #     anneal_epsilon_r = anneal_epsilon_r*e_annealing
# #     dist_to0 = min(abs(r),epsilon = anneal_epsilon_r)               # the learning rate needs to be small near 0 to avoid blowing up there.
# #     anneal_temperature = (anneal_temperature+1-temperature)*T_annealing + temperature-1
# #     r = SGLD(r, distribution_snookerscaler, dist_to0, tau = anneal_temperature)
# #     R = rbind(R,r)
# #   }
# # }
# for (t in 1:(epochR)){
#   for (k in 1:nBatch){
#     Batch = Data[which(Data$augment == k),]
#     temp <- function(theta){
#       temp = nonlinearReg(theta, Batch, nData)
#     }
#     distribution_snookerscaler <- function(r){
#       distribution_snookerscaler = snooker_scaler(r,temp,theta[,paired_Index],direction_e)
#     }
#     dist_to0 = min(abs(r),epsilon = epsilon_r)               # the learning rate needs to be small near 0 to avoid blowing up there.
#     r = SGLD(r, distribution_snookerscaler, dist_to0, tau = temperature)
#     R = rbind(R,r)
#   }
# }
# 
# # #
# hist(R[1000:2000,3])


##===================================================
## =========== !!!MAIN FUNCTION!!! ==================
##===================================================
## This is the one step snooker sampler. 

## input: 
#     x:                       current positions of the population. dim(x) = [dimension, nchain]
#     temperature_select:      temperature ladder. A sequence of length nchain. 
#     Index:                   which bunch of the population we are going to update via snooker. 

## output: 
#     those states of x that got updated.
#====================
snooker_batch <-function(theta, temperature_select, Index, tau){
  paired_Index = sample_Boltzmann(theta,temperature_select,Index)     # sample xj by boltzmann weight and record the j's.
  
  direction_e  = unit_direction(theta, Index, paired_Index)     #calculate e: the unit updating direction. 
  
  r = NA
  for (attempt in 1:3) {
    try(r <- sample_r(epsilon_r, epochR, temperature = tau[paired_Index], T_annealing = 1, e_annealing = 1, Index, paired_Index = paired_Index, direction_e = direction_e))     # sample the scaler r. After consideration I choose to initialize r at N(0,1). Might need to modify.
    # I need "try" here because sometimes the sampling got screwed up with tiny chance. 
    if(!prod(is.na(r))) break       # if non of the entries of r is still NA then break
  }
  if (prod(is.na(r))) stop("r = NA")        # if r still has NA entries then stop. 
  
  
  r_matrix = t(array(r, dim = c(length(r), dimension))) 
  snooker_update = theta[,paired_Index] + r_matrix*direction_e
  return(snooker_update) 
}



