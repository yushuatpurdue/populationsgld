### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the pSGLD to produce a chain. algorithm see paper Ding(2014 google research)


###============= initialization: call functions and initial values=======================
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/mixGauss2order.R")
distribution_sgld <- mixgaussian2    

dimension     = 2
learning_rate = 0.25
alpha         = 0.2                 # tuning parameter 1
lambda        = 0.001                 # tuning parameter 2



### =======================================================
### =============== supportive functions ==================
### =======================================================

### =========function 1: rowsum_mu===========
## In order to correct the NA gradients (give the gradient a bound). Calculate the row sum of x - mu. Idea and details refer to "./SGLD.R". 

if (exists("mu")){
  rowsum_mu <- function(x){
    return(rowSums(x - mu))
  }
}


### =========function 2: onestepSGNHT===========
## Output one step update through the SGNHT scheme. 

onesteppSGLD <- function(position, V, alpha, lambda, learning_rate = 0.1){
  ## input: 
  # position        : like position of a particle in a 2nd order Langevin Dynamic. dimension*1 array. 
  # V               : kind of like the empirical fisher information. used to construct G the preconditioning matrix. 
  # alpha           : tuning parameter 1.
  # learning_rate   : stepsize of the discretized dynamic.  
  
  ## output: 
  # G               : preconditioning matrix. 
  # position_update : 1 step updated position. 
  
  dimension = dim(position)[1]
  temp      = distribution_sgld(position)        # calculate the gradient
  grad      = temp$grad 
  hess      = temp$hess 
  
  ## the following 'if' is to bound the blowed up gradients. Idea and details refer to "./SGLD.R". 
  if (sum(is.na(grad)) > 0){            # if there exists any entry being NA in grad:
    if (exists("mu")){
      NAdirection = sign(apply(position, 2, rowsum_mu))                   #  returns the direction for correction
      blowup_entry = NAdirection*(is.na(grad))                            # (1-is.na(grad)) returns which entries blowed up. sign(rowSums(x - mu)) 
      
      grad[is.na(grad)] <- 0
      grad = 200*blowup_entry + grad
      hess[is.na(hess)] <- 0
      hess = -8000*is.na(hess) + hess
      
      # print("gradient corrected")
    }
    else {
      stop("gradient blows up")
    }
  }
  
  
  # ================= update ==================
  ## 1. update V. 
  V_update = alpha*V + (1-alpha)*grad^2
  
  ## 2. find G (the proconditioning matrix) according to V_update
  G        = 1/(lambda + sqrt(V_update))
  G        = as.vector(G)
  G        = diag(G)
  
  ## 3. find the Gamma (corrector for the drift term)
  sqrtV = sqrt(V_update)
  Gamma = - (1-alpha)*grad*hess / (sqrtV * (lambda+sqrtV)^2)
  
  ## 4. calculate the drift coefficient
  drift = - G%*%(grad) + Gamma
  
  ## 5. generate the diffusion coefficient
  
  wn_length       = length(position)
  wn              = rnorm(wn_length)
  wn              = array(wn, dim = dim(position))            # generate a gaussian noise
  diffusion       = sqrt(learning_rate*G)
  
  ## 6. update position
  position_update = position + drift*learning_rate/2 + diffusion%*%wn
  
  return(list("position_update" = position_update, "V_update" = V_update))
}

### ======================================================
### =============== main function=========================
### ======================================================
pSGLD <-function(lambda, alpha, niteration, learning_rate){
  
  # position, dimension*1:
  x = runif(dimension, min = 0, max = 1)
  x = array(x, dim = c(dimension,1))
  # V, dimension*1: 
  V = 0
  V = array(V, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,niteration))   
  output[,1] = x      
  
  # run SGHMC chain
  for (t in 2:niteration) {
    temp = onesteppSGLD(x, V, alpha, lambda, learning_rate)
    x = temp$position_update
    V = temp$V_update
    
    output[,t] =x
  }
  return(output)
}

###=========================================
###=============== run function ============
###=========================================

output_pSGLD = pSGLD(lambda, alpha, 200000, learning_rate)

plot(output_pSGLD[1,1:50000], output_pSGLD[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "pSGLD: 50000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

plot(output_pSGLD[1,], output_pSGLD[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "pSGLD: 200000 iteration", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

write.csv(output_pSGLD, "./variousSGLD/output_pSGLD.csv")



##====================================
## =========== test ==================
##====================================

### test for onesteppSGLD applied on 
V = array(0, dim = c(dimension,1))




