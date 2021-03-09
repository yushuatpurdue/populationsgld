### This file is created by Yushu Huang for EMC-SGLD experiments.
### This is the one step MSGLD sampler. See paper Kim(2020)
### An additional supportive function is defined to deal with the gradient blowup. 

###============= initialization: call functions and initial values=======================
library(MASS)
setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian   

dimension = 2
learning_rate = 0.021
temperature = 2
smoothingfactor = 0.5                        # smoothing factor beta (tailing off coefficient) for the momentum wrt time. 
biasfactor = 0.5                            # bias factor: how much portion of momentum will be used in updating position. 



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


### =========function 2: onestepMSGLD===========
## Output one step update through the SGNHT scheme. 

onestepMSGLD <- function(position, momentum, smoothingfactor, biasfactor, learning_rate = 0.1, temperature){
  ## input: 
  # position       : like position of a particle in a 2nd order Langevin Dynamic. dimension*1 array. 
  # momentum       : like momentum of a particle in a 2nd order Langevin Dynamic. dimension*1 array. 
  # friction       : kind of like the friction of the heatbath (not exactly). However it's changing wrt momentum. scaler. 
  # learning_rate  : stepsize of the discretized dynamic.  
  
  dimension = dim(position)[1]
  grad = distribution_sgld(position)$grad        # calculate the gradient
  
  ## the following 'if' is to bound the blowed up gradients. Idea and details refer to "./SGLD.R". 
  if (sum(is.na(grad)) > 0){            # if there exists any entry being NA in grad:
    if (exists("mu")){
      NAdirection = sign(apply(position, 2, rowsum_mu))                   #  returns the direction for correction
      blowup_entry = NAdirection*(is.na(grad))         # (1-is.na(grad)) returns which entries blowed up. sign(rowSums(x - mu)) 
      
      grad[is.na(grad)] <- 0
      grad = 200*blowup_entry + grad
      
      # print("gradient corrected")
    }
    else {
      stop("gradient blows up")
    }
  }
  
  wn_length = length(position)
  wn = rnorm(wn_length)
  wn = array(wn, dim = dim(position))            # generate a gaussian noise
  
  # ================= update ==================
  position_update = position - learning_rate*(grad + biasfactor*momentum) + sqrt(2*temperature*learning_rate)*wn
  momentum_update = smoothingfactor*momentum + (1 - smoothingfactor)*grad
  
  
  return(list("position_update" = position_update, "momentum_update" = momentum_update))
}

### ======================================================
### =============== main function=========================
### ======================================================

MSGLD <- function(smoothingfactor, biasfactor, temperature, niteration, learning_rate){
  # position, dimension*1:
  x = runif(dimension, min = 0, max = 1)
  x = array(x, dim = c(dimension,1))
  # momentum, dimension*1: 
  m = 0
  m = array(m, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,niteration))   
  output[,1] = x      
  
  # run SGHMC chain
  for (t in 2:niteration) {
    
    temp       = onestepMSGLD(x, m, smoothingfactor, biasfactor, learning_rate, temperature)
    x          = temp$position_update
    m          = temp$momentum_update
    output[,t] = x
    
    # print(x)
    # print(t)
  }
  return(output)
}


###=========================================
###=============== run function ============
###=========================================

output_MSGLD = MSGLD(smoothingfactor, biasfactor, temperature, 200000, learning_rate)

plot(output_MSGLD[1,1:50000], output_MSGLD[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "MSGLD: 50000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

plot(output_MSGLD[1,], output_MSGLD[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "MSGLD: 200000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

write.csv(output_MSGLD, "./variousSGLD/output_MSGLD.csv")


##====================================
## =========== test ==================
##====================================

### test for SGNHT applied on 

# position = mvrnorm(1, mu = rep(0,dimension), Sigma = diag(1,dimension))
# position = array(position,dim = c(2,1))
# momentum = array(1,dim = c(2,1))
# friction = 1
## x = array(c(1,-0.5,0.3,0.2,0.4),dim = c(1,5))
# learning_rate = 0.1
# onestepSGNHT(position, momentum, friction, A, learning_rate = 0.1)

