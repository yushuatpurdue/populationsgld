### This file is created by Yushu Huang for local SPOS experiments.
### This file runs SPOS with localized K. Algorithm improves Zhang 2020 SPOS paper. 

library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

# set.seed(1)

dataset = 'datasetNLR3'
Data = read.csv(paste0("dataset/",dataset,".csv"))
print(dataset)
source("potential/nonlinearReg.R")
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)

## =============== function: calculate K matrix ====================
# K is an M*M matrix, M = # of particles.
# K = [[k11,k12,k13,...,k1M],[k21,k22,k23,...,k2M],...,[k1M,k2M,k3M,...kMM]]. Note: row and column is flipped against convention. 

kernelComp <- function(x,eta){
  ## Input: 
  #     x       particles.      d*M matrix. d = dimension, M = number of particles.
  #     eta     scaling parameter. scale.
  ## output:
  #     K       kernel matrix. Matrix of dimension M*M. 
  #     gradK   gradient (of second term) of K matrix, summed over the second tem. 
  # Detail about gradK: it is the summation of the 4th term in Zhang 2020 equation (9). 
  # K is the kernel matrix where Kij = exp(-(xi - xj)^2/eta^2). gradK first takes gradient K wrt. xj, and this provides an M*M*d tensor. Then it sums out wrt j, and this provides an M*d matrix. 
  ## ====== compute K ===========
  M = dim(x)[2]     # M is the number of particles.
  K = array(NA, dim = c(M,M)) 
  ones = t(rep(1,M))
  ## idea: for every particle, copy their position as a matrix and take difference to the rest of the particles. only do the upper triangle to save half computation. 
  for (i in 1:M){
    xi = x[,i]%*% ones
    diffx = xi - x[,i:M]
    l2 = colSums(diffx^2)
    result = exp(-l2/eta^2)
    K[i:M,i] = result
    K[i,i:M] = result
    ones = t(ones[-1])
  }
  
  ## ======= compute gradK =========
  gradK = array(NA, dim = c(dimension,M)) 
  ones = t(rep(1,M))
  ## idea: for every particle, copy their position as a matrix and take difference to the rest of the particles. only do the upper triangle to save half computation. 
  for (i in 1:M){
    xi = x[,i]%*% ones
    diffx = xi - x
    l2 = 2*diffx/eta^2
    result = l2 %*% K[,i]
    gradK[,i] = rowSums(result)
  }
  return(list('K' = K,'gradK' = gradK))
}

### =========function: rowsum_mu===========
## In order to correct the NA gradients (give the gradient a bound). Calculate the row sum of x - mu. Idea and details refer to "./SGLD.R". 

if (exists("mu")){
  rowsum_mu <- function(x){
    return(rowSums(x - mu))
  }
}

## =============== function: one step update SPOS ====================

onestepSPOS <- function(x, lr, eta, Beta){
  ## Input:
  #   x    positions of the particles. d*M, d dimension, M #of particles.
  #   lr   learning rate
  #   eta  scaling parameter in kernel
  #   beta tuning parameter. Represents temperature inversed. 
  ## Side note: Can I set different temperature for each particle? (In paper it is the same for all. )
  ## Output:
  #   one step update of the positions. 
  
  M = dim(x)[2]     # M is the number of particles.
  ## calculate gradient:
  gradX = potential(x)$grad
  
  ## the following 'if' is to bound the blowed up gradients. Idea and details refer to "./SGLD.R". 
  if (sum(is.na(gradX)) > 0){            # if there exists any entry being NA in grad:
    if (exists("mu")){
      NAdirection = sign(apply(x, 2, rowsum_mu))                   #  returns the direction for correction
      blowup_entry = NAdirection*(is.na(gradX))         # (1-is.na(grad)) returns which entries blowed up. sign(rowSums(x - mu)) 
      
      gradX[is.na(gradX)] <- 0
      gradX = 200*blowup_entry + gradX
      
      # print("gradient corrected")
    }
    else {
      stop("gradient blows up")
    }
  }
  
  ## SGLD term:
  wn_length = length(x)
  wn        = rnorm(wn_length)
  wn        = array(wn, dim = dim(x))            # generate a gaussian noise
  sgldUpdate = x - lr*gradX/Beta + sqrt(2*lr/Beta)*wn
  
  ## variational gradient
  K_gradK = kernelComp(x,eta)
  K = K_gradK$K
  gradK = K_gradK$gradK
  VG1 = -gradX%*%t(K)*lr/M + gradK*lr/M
  
  xUpdate = sgldUpdate + VG1
  return(xUpdate)
}


SPOSbatch <- function(nepoch, nchain, learning_rate, eta, Beta){
  ###=========== set path ================
  # if(getwd() != "/Users/yushuhuang/Documents/pu/2020fall/populationsgld"){
  #   setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  # }
  
  ###============= initialization: call functions and initial values=======================
  
  theta <<- t(mvrnorm(nchain, mu = c(0,0), Sigma = diag(1,2)))
  
  ###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,nepoch*nBatch))   # result array
  # output[,,1] = x                               # record initialization
  
  for (t in 1:nepoch) {
    for (k in 1:nBatch) {
      Batch = Data[which(Data$augment == k),]
      potential <<- function(theta){
        potential = stochPotential(theta, Batch, nData)
      }
      theta          = onestepSPOS(theta, learning_rate, eta, Beta)
      
      output[,,((t-1)*nBatch+k)] <- theta
      # print(theta)
    }
    
    
    
    # print(t)
  }
  
  return(list("output" = output))
  
}


# ============ test =================

output = SPOSbatch(nepoch = 40, nchain=100, learning_rate=0.03, eta=2, Beta=5)
rowMeans(output$output[,1,3000:4000])
plot(output$output[1,1,3000:4000], output$output[2,1,3000:4000],pch = ".")

