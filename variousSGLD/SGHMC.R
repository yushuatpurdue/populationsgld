### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the SGHMC to produce a chain. See paper Chen (2014 ICML). 
### An additional supportive function is defined to deal with the gradient blowup. 

###============= initialization: call functions and initial values=======================
library(MASS)
setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian    

## initialize parameters: 
dimension = 2
M         = 1                  # mass of the particle. can be a dimension*dimension matrix, but need to modify the code in "onestepSGHMC()". 
C         = 14                # tuning parameter, friction. 
learning_rate = 0.1
temperature = 5.6


### ======================================================
### =============== supportive functions==================
### ======================================================

### =========function 1: : rowsum_mu===========
## In order to correct the NA gradients (give the gradient a bound). Calculate the row sum of x - mu. Idea and details refer to "./SGLD.R". 

if (exists("mu")){
  rowsum_mu <- function(x){
    return(rowSums(x - mu))
  }
}


### =========function 2: onestepSGHMC===========
## Output one step update through the SGHMC scheme. 

onestepSGHMC <- function(position, momentum, learning_rate = 0.1, temperature = 1 ){
  
  position_update = position + momentum/M*learning_rate
  grad = distribution_sgld(position_update)$grad      # compute the gradient. 
  
  ## the following 'if' is to bound the blowed up gradients. Idea and details refer to "./SGLD.R". 
  if (sum(is.na(grad)) > 0){            # if there exists any entry being NA in grad:
    if (exists("mu")){
      NAdirection = sign(apply(position_update, 2, rowsum_mu))                   #  returns the direction for correction
      blowup_entry = NAdirection*(is.na(grad))         # (1-is.na(grad)) returns which entries blowed up. sign(rowSums(x - mu)) 
      
      grad[is.na(grad)] <- 0
      grad = 200*blowup_entry + grad
      
      # print("gradient corrected")
    }
    else {
      stop("gradient blows up")
    }
  }
  
  wn_length = length(momentum)
  wn = rnorm(wn_length)
  wn = array(wn, dim = dim(momentum))

  momentum_update = momentum - t(t(grad) * learning_rate) - momentum*C/M*learning_rate + sqrt(2*C*learning_rate*temperature)*wn
  
  return(list("position_update" = position_update, "momentum_update" = momentum_update))
}



### ======================================================
### =============== main function=========================
### ======================================================

SGHMC <- function(temperature, niteration, learning_rate){
  # position:
  x = runif(dimension, min = 0, max = 1)
  x = array(x, dim = c(dimension,1))
  # momentum: 
  m = 0
  m = array(m, dim = c(dimension,1))
  
  # result array:
  output = array(NA,dim = c(dimension,2,niteration))   
  output[,1,1] = x      
  output[,2,1] = m
  
  # run SGHMC chain
  for (t in 2:niteration) {
    
    temp = onestepSGHMC(x, m, learning_rate, temperature)
    x = temp$position_update
    m = temp$momentum_update
    output[,1,t] = x
    output[,2,t] = m
    
    # print(x)
    # print(t)
  }
  return(output)
}



###=========================================
###=============== run function ============
###=========================================

output_SGHMC = SGHMC(temperature, 200000, learning_rate)

plot(output_SGHMC[1,1,1:50000], output_SGHMC[2,1,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGHMC: 50000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

plot(output_SGHMC[1,1,], output_SGHMC[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGHMC: 200000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

write.csv(output_SGHMC, "./variousSGLD/output_SGHMC.csv")

##====================================
## =========== test ==================
##====================================

### test for SGHMC applied on 

# position = mvrnorm(1, mu = rep(0,dimension), Sigma = diag(1,dimension))
# position = array(position,dim = c(2,1))
# momentum = array(c(0,0),dim = c(2,1))
# # x = array(c(1,-0.5,0.3,0.2,0.4),dim = c(1,5))
# learning_rate = 0.01
# # onestepSGHMC(r,epsilon)





