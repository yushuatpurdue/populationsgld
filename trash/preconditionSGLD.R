## This is the preconditioned SGLD sampler. 

## input: 
#      r0:            starting positions of the chains. 
# distribution_sgld:  target distribution. 
# niteration:         number of iterations
# epsilon:            learning rate. default 0.1. (SGLD)
#    tau:             temperature. default all 1. (SGLD)
# lambda:             tunning parameter, default 1e-5. (precondition SGLD)
# alpha:              tuning parameter, default 0.99.  (precondition SGLD)

## output: 
#      A chain of positions. 

###=========== redundant block: should have realized in another file where this file is called ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/2gauss_hess.R")                   # source some potential
# distribution_sgld <-  twogaussian                   # pass it to a general name that is recognizable for function SGLD. 
# rm(twogaussian)                                     # remove the redundant source function. 

### ====================================================
### ================ draft ===================
### ====================================================

## ======================== one step pSGLD =============================

## initialization: 
V = rep(0,5)
r0 = runif(5, -2,2)
pd_grad_hess_r = snooker_scaler(r, twogaussian_hess,paired_x, direction)
grad_r = pd_grad_hess_r$grad 
hess_r = pd_grad_hess_r$hess 

## from current V and r to the next V and r. 
V = alpha*V + (1-alpha)*grad_r^2
sqrtV = sqrt(V)
G = 1/(lambda + sqrtV)
Gamma_r = -(1 - alpha)*G^2/sqrtV*grad_r*hess_r

r = r - epsilon/2*(G*grad_r - Gamma_r)
r = r + mvrnorm(1, rep(0,nsample), diag(epsilon*G))




### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================

preconditionSGLD <- function(r0, niteration, distribution_sgld, epsilon = 0.1,tau = 1, alpha = 1e-5, lambda = 0.99){
  
  nsample = length(r0)
  V = rep(0,nsample)
  r_chain = rbind(r0)
  
  for (i in 1:niteration-1){
    r = r_chain[nrow(r_chain),]
    pd_grad_hess_r = distribution_sgld(r, twogaussian_hess,paired_x, direction)
    grad_r = pd_grad_hess_r$grad 
    hess_r = pd_grad_hess_r$hess 
    
    ## from current V and r to the next V and r. 
    V = alpha*V + (1-alpha)*grad_r^2
    sqrtV = sqrt(V)
    G = 1/(lambda + sqrtV)
    Gamma_r = -(1 - alpha)*G^2/sqrtV*grad_r*hess_r
    
    r = r - epsilon/2*(G*grad_r - Gamma_r)
    r = r + mvrnorm(1, rep(0,nsample), diag(epsilon*G))
    
    r_chain = rbind(r_chain,r)
    
  }
  
  return(r_chain)
}



##====================================
## =========== test ==================
##====================================

### 
source("potential/snooker_scaler.R")
source("potential/2gauss_hess.R")
r0 = runif(5, -2,2)
# paired_x = t(mvrnorm(5, mu = rep(0,2), Sigma = diag(1,2)))
paired_x = x[,paired_Index]
direction = direction_e
result = preconditionSGLD(r0, 10000, snooker_scaler)
hist(result[,1], breaks = seq(-12,12,0.2))

