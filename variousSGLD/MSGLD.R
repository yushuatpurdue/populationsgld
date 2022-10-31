### This file is created by Yushu Huang for EMC-SGLD experiments.
### This is the one step MSGLD sampler. See paper Kim(2020)
### An additional supportive function is defined to deal with the gradient blowup. 

###============= initialization: call functions and initial values=======================
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/mixtureGaussian.R")
# distribution_sgld <- mixgaussian   
# 
dimension = 2
learning_rate = 0.001
temperature = 1
smoothingfactor = 0.9                        # smoothing factor beta (tailing off coefficient) for the momentum wrt time. beta_1 in pater
biasfactor = 28                            # bias factor: how much portion of momentum will be used in updating position. alpha in paper

## ========= below is the initialization for minibatch sampling narrow potential ravine example ======

Data = read.csv("dataset/datasetNLR1.csv")
source("potential/nonlinearReg.R")
# potential = nonlinearReg
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)




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
  # smoothingfactor: how much past momentum is cumulated into current momentum.
  # biasfactor     : how much momentum is used for updating. 
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

### =========function 1: MSGLD===========
MSGLD <- function(smoothingfactor, biasfactor, temperature, niteration, learning_rate){
  ## input: 
  # smoothingfactor: how much past momentum is cumulated into current momentum.
  # biasfactor     : how much momentum is used for updating. 
  # temperature    : temperature of langevin dynamic
  # niteration     : number of epoch
  # learning_rate  : stepsize of the discretized dynamic. 
  
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

### =========function 2: batchMSGLD===========

batchMSGLD <- function(smoothingfactor, biasfactor, temperature, nepoch, learning_rate){
  ## input: 
  # smoothingfactor: how much past momentum is cumulated into current momentum.
  # biasfactor     : how much momentum is used for updating. 
  # temperature    : temperature of langevin dynamic
  # nepoch         : number of iteration
  # learning_rate  : stepsize of the discretized dynamic.  
  
  # position, dimension*1:
  # x = runif(dimension, min = 0, max = 1)
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # momentum, dimension*1: 
  m = 0
  m = array(m, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,nepoch*nBatch))   
  output[,1] = x      
  
  # run SGHMC chain
  
  for (t in 1:nepoch) {
    for (k in 1:nBatch) {
      Batch = Data[which(Data$augment == k),]
      distribution_sgld <<- function(theta){
        distribution_sgld = stochPotential(theta, Batch, nData)
      }
      temp       = onestepMSGLD(x, m, smoothingfactor, biasfactor, learning_rate, temperature)
      x          = temp$position_update
      m          = temp$momentum_update
      output[,((t-1)*nBatch+k)] = x
    }
  }
  
  
  return(output)
}



###=========================================
###=============== run function ============
###=========================================

output_MSGLD = batchMSGLD(smoothingfactor, biasfactor, temperature, 2000, learning_rate)
# 
# plot(output_MSGLD[1,1:50000], output_MSGLD[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "MSGLD: 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output_MSGLD[1,], output_MSGLD[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "MSGLD: 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")

# write.csv(output_MSGLD, "./variousSGLD/output_MSGLD.csv")

# nrun    = 10                  # run nrun independent times.
# mu_est_MSGLD = array(NA, dim = c(2, nrun))
# var_est_MSGLD = array(NA, dim = c(2,2,nrun))
# 
# for (i in 1:nrun) {
#   print(paste0(i,"th run"))
#   start_time <- Sys.time()
#   output_MSGLD = MSGLD(smoothingfactor, biasfactor, temperature, 150000000, learning_rate)
#   end_time <- Sys.time()
#   print(paste0("running time = ", end_time - start_time))
#   # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
#   # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))
# 
# 
#   record = output_MSGLD[,c(1:50000)*3000]
# 
#   mu_est_MSGLD[,i] = rowMeans(record)
#   var_est_MSGLD[,,i] = cov(t(record))
# 
# 
#   pdf(file=paste0("mixGaussianMSGLDplot_",i,".pdf"))
#   for (j in 1:50) {
#     plot(record[1,c(1:(1000*j))], record[2,c(1:(1000*j))], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("MSGLD ",1000*j*3000, " iteration"), xlab = "x", ylab = "y")
#     points(mu[1,], mu[2,], pch = 3,col = "green")
# 
#   }
#   plot(record[1,],record[2,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "MSGLD sample path", xlab = "x", ylab = "y")
# 
#   dev.off()
#   print(paste("mu estimate = ", mu_est_MSGLD[,i]))
#   print(paste("var estimate = ", var_est_MSGLD[,,i]))
# }
# 
# mu_est = rowMeans(mu_est_MSGLD)
# sigma_est_1 = mean(var_est_MSGLD[1,1,])
# sigma_est_2 = mean(var_est_MSGLD[2,2,])
# cov_est = mean(var_est_MSGLD[1,2,])
# print(paste0("estimated mu = ", mu_est))
# print(paste0("estimated cov = ", sigma_est_1,", ", sigma_est_2,", ",cov_est))
# 
# df = data.frame(mu1 = mu_est_MSGLD[1,], mu2 = mu_est_MSGLD[2,], var1 = var_est_MSGLD[1,1,], var2 = var_est_MSGLD[2,2,], cov = var_est_MSGLD[1,2,]) # estimating mu and sigma.
# write.csv(df, "./mu_var_est_MSGLD.csv")


