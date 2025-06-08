### This file is created by Yushu Huang for local SPOS experiments.
### This file runs SPOS with localized K. Algorithm improves Zhang 2020 SPOS paper. 
### It also adds a temperature ladder. 

library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

set.seed(1)

dataset = 'datasetNLR1'
Data = read.csv(paste0("SPOS/dataset/",dataset,".csv"))
print(dataset)
source("potential/nonlinearReg.R")
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)


## =============== function: calculate K matrix ====================
# K is an M*M matrix, M = # of particles.
# K = [[k11,k12,k13,...,k1M],[k21,k22,k23,...,k2M],...,[k1M,k2M,k3M,...kMM]]. Note: row and column is flipped against convention. 

truncKernelComp <- function(x,eta,truncAt){
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
  K[K<truncAt] <- 0
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

onestepLSPOS <- function(x, lr, eta, Beta, temperature, truncAt){
  ## Input:
  #   x    positions of the particles. d*M, d dimension, M #of particles.
  #   lr   learning rate
  #   eta  scaling parameter in kernel
  #   beta tuning parameter. Represents temperature. 
  ## Side note: Can I set different temperature for each particle? (In paper it is the same for all. )
  ## Output:
  #   one step update of the positions. 
  
  M = dim(x)[2]     # M is the number of particles.
  ## calculate gradient:
  gradX = potential(x)$grad
  lr = t(array(lr,dim = dim(t(x))))
  Beta = t(array(Beta,dim = dim(t(x))))
  temperature = t(array(temperature, dim = dim(t(x))))
  
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
  sgldUpdate = x - lr*gradX/Beta + sqrt(2*lr*temperature/Beta)*wn
  
  ## variational gradient
  K_gradK = truncKernelComp(x,eta,truncAt)
  K = K_gradK$K
  gradK = K_gradK$gradK
  VG1 = -gradX%*%t(K)*lr/M + gradK*lr/M
  
  xUpdate = sgldUpdate + VG1
  return(xUpdate)
}


LSPOSbatch <- function(nepoch, nchain, learning_rate, eta, Beta, temperature, truncAt){
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
      theta          = onestepLSPOS(theta, learning_rate, eta, Beta, temperature, truncAt)
      
      output[,,((t-1)*nBatch+k)] <- theta
      # print(theta)
    }
    
    
    
    # print(t)
  }
  
  return(list("output" = output))
  
}



# ============ test =================

nrun = 40
M = 500
# lr = seq(0.03,0.06,0.03/(M - 1))
# temperature = seq(1,20,2/(M - 1))
eta = 3
Beta = 50
lr = 0.5
temperature = 1
cutoff = exp(-0.25)
# cutoff = 0   # SPOS

mu_est = array(NA, dim = c(2, nrun))

for (seed in c(1:nrun)){
  print(paste0("LSPOS seed = ", seed))
  set.seed(seed)
  output = LSPOSbatch(nepoch = 40, nchain=M, learning_rate=lr, eta=eta, Beta=Beta, temperature= temperature, truncAt = cutoff)
  pdf(file=paste0("./SPOS/LSPOS",seed,".pdf"))
  plot(output$output[1,1,1000:4000], output$output[2,1,1000:4000],pch = ".", main = paste0("eta=",eta,",beta=", Beta,",lr=",lr, ",cutoff=", cutoff))
  dev.off()
  mean_est = rowMeans(output$output[,1,1000:4000])
  mu_est[,seed] = mean_est
  print(paste0("LSPOS mean est = ", mean_est))
}
write.csv(mu_est, paste0("./SPOS/thetaEstLSPOS_M",M,".csv"))

# SPOS 
cutoff = 0   # SPOS
mu_est = array(NA, dim = c(2, nrun))
for (seed in c(1:nrun)){
  print(paste0("SPOS seed = ", seed))
  set.seed(seed)
  output = LSPOSbatch(nepoch = 40, nchain=M, learning_rate=lr, eta=eta, Beta=Beta, temperature= temperature, truncAt = cutoff)
  pdf(file=paste0("./SPOS/SPOS",seed,".pdf"))
  plot(output$output[1,1,1000:4000], output$output[2,1,1000:4000],pch = ".", main = paste0("eta=",eta,",beta=", Beta,",lr=",lr, ",cutoff=", cutoff))
  dev.off()
  mean_est = rowMeans(output$output[,1,1000:4000])
  mu_est[,seed] = mean_est
  print(paste0("SPOS mean est = ", mean_est))
}
write.csv(mu_est, paste0("./SPOS/thetaEstSPOS_M",M,".csv"))

# set.seed(9)
# output = LSPOSbatch(nepoch = 40, nchain=M, learning_rate=lr, eta=eta, Beta=Beta, temperature= temperature, truncAt = cutoff)
# plot(output$output[1,1,1000:4000], output$output[2,1,1000:4000],pch = ".", main = paste0("eta=",eta,",beta=", Beta,",lr=",lr, ",cutoff=", cutoff))
# rowMeans(output$output[,1,1000:4000])


thetaEstLSPOS_M500 = read.csv("/Users/yushuhuang/Documents/pu/research/populationsgld/SPOS/thetaEstLSPOS_M500.csv")
thetaEstSPOS_M500 = read.csv("/Users/yushuhuang/Documents/pu/research/populationsgld/SPOS/thetaEstSPOS_M500.csv")
