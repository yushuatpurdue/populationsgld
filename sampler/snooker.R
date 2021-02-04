### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the one step snooker sampler. 

## input: 
#     x:       current positions of the population. dim(x) = [dimension, nsample]
#     Index:   which bunch of the population we are going to update via snooker. 

## output: 
#     those states of x that got updated.



##==============================================
## =========== initialization ==================
##==============================================


###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

### ========== first: parameter setups ===================================================
nchain = dim(x)[2]                    # number of chains
Ts      = seq(1,5,4/(nchain-1))            # selection temperature for sampling xj
### ========== call function: calculating potential & gradient ===============

source("potential/mixtureGaussian.R")  
epsilon_r = min(0.1, sigma^2)                        # learning rate for sgld in sampling r. 
potential = mixgaussian

# ==================================
## load the functions to be used when sampling. 

source("potential/snooker_scaler.R")              # source the potential of r. 

source("sampler/SGLD.R")                          # source 1 step SLGD as updating method. 



##==================================================
## =========== necessary function ==================
##==================================================

### ============ Second: find the chains to pair with x ===========
# ==================================
## function calculate the boltzmann weight. 

Boltzmann <- function(x, Ts){
  # input:  x            = vector of positions (again vector), so overall matrix
  #         Ts           = temperature selected for sampling, scaler
  # output: Boltzmann(x) = their boltzmann weight = exp(-u(x)/Ts), vector
  pd_grad = potential(x)
  Boltzmann = pd_grad$pd / Ts
  return(Boltzmann)
}

# ==================================
## sample xj by boltzmann weight. 
sample_Boltzmann <- function(x, Ts, index){
  # input:  x            = vector of positions (again vector), so overall matrix
  #         Ts           = temperature selected for sampling, scaler
  #         index        = the indices that should be take away before sampling
  # output: sampleid     = the indices of the sampled chain in {x} to snooker with xi
  sample_Boltzmann = c()
  for (i in index) 
  {
    weight = Boltzmann(x, Ts)     # calculate the boltzmann weights for everyone in {x}
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
  diff_pair = x[,index] - x[,pair]
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
    r = SGLD(r, distribution_snookerscaler, dist_to0)
    # R = rbind(R,r)
  }
  return(r)
}

## test for the convergence of sampling r.
#   # r0 = runif(5, -1,1)
#   # r0[which(r0<0)] = r0[which(r0<0)]*5-0.5
#   # r0[which(r0>0)] = r0[which(r0>0)]*5+0.5
#   r0 = rnorm(length(Index))
#   r0 = array(r0,dim = c(1,length(r0)))
#   R = matrix(ncol = 5)
#   R[1,] = r0
#   r = r0
#   for (i in 1:5000) {
#     dist_to0 = min(abs(r),epsilon_r)
#     r = SGLD(r,distribution_snookerscaler, dist_to0)
#     R = rbind(R,r)
#   }
# # #
# hist(R[,5], breaks = seq(-12,12,0.2))


##===================================================
## =========== !!!MAIN FUNCTION!!! ==================
##===================================================

#====================
snooker <-function(x, Ts, Index){
  paired_Index = sample_Boltzmann(x,Ts,Index)     # sample xj by boltzmann weight and record the j's.
  direction_e  = unit_direction(x, Index, paired_Index)     #calculate e: the unit updating direction. 
  r = rnorm(length(Index))
  r = array(r,dim = c(1,length(r)))
  
  # source("potential/2gauss_hess.R")
  
  distribution_snookerscaler <- function(r){
    distribution_snookerscaler = snooker_scaler(r,potential,x[,paired_Index],direction_e)
  }
  
  r = sample_r(epsilon_r, Index, distribution_snookerscaler)     # sample the scaler r. After consideration I choose to initialize r at N(0,1). Might need to modify. 

  r_matrix = t(array(r, dim = c(length(r), dimension))) 
  snooker_update = x[,paired_Index] + r_matrix*direction_e 
  return(snooker_update) 
}



##====================================
## =========== test ==================
##====================================
# nchain = 10
# dimension = 2
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# #
# Index = sample.int(10,5)        # sample 5 integers between 1 and 10 without replacement.
# paired_Index = sample.int(9,5)  # sample j so xj that gotta paired with xi (just an example. practically need to use boltzmann weight)
# # if j>i, means the real j is j+1 (because j is in 1:9 and need to be in 1:10 while jump over i). 
# addon = rep(0,nchains)        
# addon = addon(which(paired_Index>Index))
# paired_Index = paired_Index + addon
# 
# for (i in Index) {
#   print(i)
# }

