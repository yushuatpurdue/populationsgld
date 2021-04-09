### This file is created by Yushu Huang for population SGMCMC experiments.
##  This is the one step snooker sampler. Main function named snooker, and a bunch of supportive functions. 
##  Follows the EMC algorithm (Faming Liang 2001) and reconstruct it under the SGLD setting. 

##=======================================================
## =========== structural introduction ==================
##=======================================================
## THIS CODE CONSISITS OF 3 PARTS. 

## part 1: initialization

## part 2: supportive functions. 
  # Boltzmann        : calculating the boltzmann weights of the current chains, for the subsequent sampling. 
  # sample_Boltzmann : sample the paired chains xj according to the Boltzmann weight. 
  # l2               : calculate the l2 distance, to normalize xj - xi
  # unit_direction   : calculate normalized direction of xj - xi
  # sample_r         : sample the snooker scaler. 

## part 3: main function: snooker
  # for a population of chains, uniformly draw a chain xi (this step is not included in this file), and pair it with xj chain drawn under Boltzmann distribution. 
  # calculate the unit direction, and sample the snooker scaler. 
  # updated xi = xj + r*e 




##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: parameter setups ===================================================
nchain  = dim(x)[2]                    # number of chains
# temperature_select = 0.5             # selection temperature for sampling xj
# Ts      = seq(1,5,4/(nchain-1))      # temperature ladder. Actually no use in snooker...
### ========== call function: calculating potential & gradient ===============

source("potential/mixtureGaussian.R")                # source the potential & gradient calculator. 
epsilon_r = min(0.1, 2*sigma^2)                        # learning rate for sgld in sampling r. 
potential = mixgaussian

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
sample_r <- function(epsilon,Index, distribution_snookerscaler){
  r = rnorm(length(Index))
  r = array(r,dim = c(1,length(r)))
  # R = matrix(ncol = 5)
  # R[1,] = r
  for (i in 1:5000) {
    dist_to0 = min(abs(r),epsilon = epsilon_r)               # the learning rate needs to be small near 0 to avoid blowing up there.  
    r = SGLD(r, distribution_snookerscaler, dist_to0, tau = 2)
    # R = rbind(R,r)
  }
  return(r)
}

# test for the convergence of sampling r.
  # r0 = runif(10, -1,1)
  # r0[which(r0<0)] = r0[which(r0<0)]*5-0.5
  # r0[which(r0>0)] = r0[which(r0>0)]*5+0.5
#   r0 = rnorm(length(Index))
#   r0 = array(r0,dim = c(1,length(r0)))
#   R = matrix(ncol = 10)
#   R[1,] = r0
#   r = r0
#   for (i in 1:1000) {
#     dist_to0 = min(abs(r),0.02)
#     r = SGLD(r,distribution_snookerscaler, dist_to0, tau = 2)
#     R = rbind(R,r)
#   }
# # #
# hist(R[,7], breaks = seq(-12,12,0.2))


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
snooker <-function(x, temperature_select, Index){
  paired_Index = sample_Boltzmann(x,temperature_select,Index)     # sample xj by boltzmann weight and record the j's.
  direction_e  = unit_direction(x, Index, paired_Index)     #calculate e: the unit updating direction. 
  
  # source("potential/2gauss_hess.R")
  
  distribution_snookerscaler <- function(r){
    distribution_snookerscaler = snooker_scaler(r,potential,x[,paired_Index],direction_e)
  }
  
  r = NA
  for (attempt in 1:3) {
    try(r <- sample_r(epsilon_r, Index, distribution_snookerscaler), silent = TRUE)     # sample the scaler r. After consideration I choose to initialize r at N(0,1). Might need to modify.
    # I need "try" here because sometimes the sampling got screwed up with tiny chance. 
    if(!prod(is.na(r))) break       # if non of the entries of r is still NA then break
  }
  if (prod(is.na(r))) stop("r = NA")        # if r still has NA entries then stop. 
  
  
  r_matrix = t(array(r, dim = c(length(r), dimension))) 
  snooker_update = x[,paired_Index] + r_matrix*direction_e 
  return(snooker_update) 
}



##====================================
## =========== test ==================
##====================================
# nchain = 20
# dimension = 2
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# #
# Index = sample.int(20,10)        # sample 10 integers between 1 and 20 without replacement.
# paired_Index = sample.int(9,5)  # sample j so xj that gotta paired with xi (just an example. practically need to use boltzmann weight)
# # if j>i, means the real j is j+1 (because j is in 1:9 and need to be in 1:10 while jump over i). 
# addon = rep(0,nchains)        
# addon = addon(which(paired_Index>Index))
# paired_Index = paired_Index + addon
# 
# for (i in Index) {
#   print(i)
# }

