### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the pSGLD to produce a chain. algorithm see paper Ruqi Zhang(2019 ICLR)


###============= mixgaussian initialization: call functions and initial values=======================
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")

# source("potential/mixtureGaussian.R")
# distribution_sgld <- mixgaussian

set.seed(1)
source("sampler/SGLD.R")                          # source 1 step SGLD as part of updating method. 

dimension     = 2
learning_rate = 0.002
K = 30000      # number of total iterations
M = 300         # number of cycles
beta = 0.1    # proportion of exploration stage
Ts = 1         # temperature


## ========= below is the initialization for minibatch sampling narrow potential ravine example ======

Data = read.csv("dataset/datasetNLR5.csv")
source("potential/nonlinearReg.R")
# potential = nonlinearReg
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)

### =======================================================
### =============== supportive functions ==================
### =======================================================

### =========function 1: learning rate computation===========

lrCompute <-function(currentIter, K,M,alpha0){
  ## input: 
  # currentIter     : Now we compute learning rate for currentIter'th iteration.  
  # K               : number of total iterations 
  # alpha0          : initial learning rate.
  # M               : number of cycles. 
  
  ## output: 
  # lrCompute       : learning rate of iteration currentIter. 
  # propCycle       : the proportion of a cycle that the sampler has reached. 
  
  iterPerCycle = ceiling(K/M)
  propCycle = (currentIter-1)%%iterPerCycle/iterPerCycle
  result = alpha0/2*(cos(pi*propCycle)+1)
  return(list("lr" = result, "propCycle" = propCycle))
  
}


### =========function 2: onestepcSGLD===========
onestepcSGLD <- function(position, currentIter, K, M, alpha0, beta){
  ## input: 
  # position        : like position of a particle in a 2nd order Langevin Dynamic. dimension*1 array.
  # currentIter     : In this step, we are at currentIter'th iteration. scaler.
  # K               : number of total iterations. scaler
  # M               : number of cycles. scaler
  # alpha0          : initial learning rate. scaler
  # beta            : proportion in a cycle where exploration stage is at. scaler
  
  ## output: 
  # position_update : 1 step updated position. 
  temp     = lrCompute(currentIter, K,M,alpha0)
  lr        = temp$lr         # learning rate at this iteration
  propCycle = temp$propCycle       # proportion in a cycle at this iteration
  
  grad      = distribution_sgld(position)$grad     # gradient
  
  if(propCycle <beta){
    position_update = position - lr* grad
  }
  else{
    position_update = SGLD(position, distribution_sgld, epsilon = lr, tau = Ts)
  }
  return(position_update)
}

### ======================================================
### =============== main function=========================
### ======================================================

### =========function 1: cSGLD===========
cSGLD <- function(niteration, K, M, alpha0, beta){
  ## input: 
  # niteration      : number of iterations. scaler
  # K               : number of total iterations. scaler
  # M               : number of cycles. scaler
  # alpha0          : initial learning rate. scaler
  # beta            : proportion in a cycle where exploration stage is at. scaler 
  
  # position, dimension*1:
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # x = runif(dimension, min = 0, max = 10)
  # x = array(x, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,niteration))   
  output[,1] = x      
  
  # run SGHMC chain
  for (t in 2:niteration) {
    
    x       = onestepcSGLD(x, t, K, M, alpha0=learning_rate, beta)
    output[,t] = x
    
    # print(x)
    # print(t)
  }
  return(output)
  
}

### =========function 2: batchcSGLD===========
batchcSGLD <- function(nepoch, K, M, alpha0, beta){
  ## input: 
  # niteration      : number of epochs. scaler
  # K               : number of total iterations. scaler
  # M               : number of cycles. scaler
  # alpha0          : initial learning rate. scaler
  # beta            : proportion in a cycle where exploration stage is at. scaler 
  
  # position, dimension*1:
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,nepoch*nBatch))   
  output[,1] = x 
  
  # run cSGLD chain
  for (t in 1:nepoch){
    for (k in 1:nBatch){
      Batch = Data[which(Data$augment == k),]
      distribution_sgld <<- function(theta){
        distribution_sgld = stochPotential(theta, Batch, nData)
      }
      x       = onestepcSGLD(x, t, K, M, alpha0=learning_rate, beta)
      output[,((t-1)*nBatch+k)] = x
    }
  }
  return(output)
}


###=========================================
###=============== run function ============
###=========================================

##=========== mixgaussian =================

# nrun    = 20                  # run nrun independent times.
# mu_est_cSGLD = array(NA, dim = c(2, nrun))
# var_est_cSGLD = array(NA, dim = c(2,2,nrun))
# niteration = 30000*30
# 
# for (i in 1:nrun) {
#   print('cycSGLD')
#   print(paste0(i,"th run"))
#   start_time <- Sys.time()
#   output_cSGLD = cSGLD(niteration, K, M, alpha0=learning_rate, beta)
#   # output_pSGLD = batchpSGLD(lambda, alpha, 50000, learning_rate = 0.001)
#   end_time <- Sys.time()
#   print(paste0("running time = ", end_time - start_time))
#   # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
#   # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))
# 
# 
#   record = output_cSGLD[,c((10000*30):niteration)]
#   mu_est_cSGLD[,i] = rowMeans(record)
#   var_est_cSGLD[,,i] = cov(t(record))
#   
#   print(paste("mu estimate = ", mu_est_cSGLD[,i]))
#   print(paste("var estimate = ", var_est_cSGLD[,,i]))
# }
# 
# print(paste0('mu est = ',rowMeans(mu_est_cSGLD)))
# print(paste0('mu est std = ',sd(mu_est_cSGLD[1,]),', ' ,sd(mu_est_cSGLD[2,])))
# print(paste0('var est = ',mean(var_est_cSGLD[1,1,]), ', ' ,mean(var_est_cSGLD[1,2,]), ', ' ,mean(var_est_cSGLD[2,1,]), ', ' ,mean(var_est_cSGLD[2,2,])))
# print(paste0('var est std = ',sd(var_est_cSGLD[1,1,]), ', ' ,sd(var_est_cSGLD[1,2,]), ', ' ,sd(var_est_cSGLD[2,1,]), ', ' ,sd(var_est_cSGLD[2,2,])))
# 
# 
# pdf(file=paste0("mixGaussiancSGLDplot_",i,".pdf"))
# for (j in 1:6) {
#   # plot(output_snooker$output[1,1,1:(1000*j)], output_snooker$output[2,1,1:(1000*j)], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("EMC ",1000*j, " iteration"), xlab = "x", ylab = "y")
#   plot(output_cSGLD[1,c(1:(5000*j))*30], output_cSGLD[2,c(1:(5000*j))*30], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("cSGLD ",5000*j*30, " iteration"), xlab = "x", ylab = "y")
#   points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# }
# # plot(output_snooker$output[1,1,],output_snooker$output[2,1,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "EMC-SGLD sample path", xlab = "x", ylab = "y")
# plot(output_cSGLD[1,c(10000:30000)*30],output_cSGLD[2,c(10000:30000)*30], type = "l", xlim = c(0,10), ylim = c(0,10), main = "cSGLD sample path", xlab = "x", ylab = "y")
# 
# dev.off()

##============== nonlinear regression ===================

nrun    = 20                  # run nrun independent times.
mu_est_cSGLD = array(NA, dim = c(2, nrun))   # nonlinear regression
nepoch = 20000

for (i in 1:nrun) {
  print('cSGLD')
  print(paste0(i,"th run"))
  start_time <- Sys.time()
  outputcSGLD = batchcSGLD(nepoch, K, M, alpha0=learning_rate, beta)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))
  print(paste0(rowMeans(outputcSGLD[,1000000:2000000])))
  
  mu_est_cSGLD[,i] = rowMeans(outputcSGLD[,1000000:2000000])
  
}
print(paste("avg mu estimate = ", rowMeans(mu_est_cSGLD)))
print(paste("mu estimator sd = ", sd(mu_est_cSGLD[1,]),', ' ,sd(mu_est_cSGLD[2,])))
