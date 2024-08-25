### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the SGHMC to produce a chain. See paper Chen (2014 ICML). 
### An additional supportive function is defined to deal with the gradient blowup. 

###============= initialization: call functions and initial values=======================
set.seed(1)
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian

## initialize parameters: 
dimension = 2
M         = 1                  # mass of the particle. can be a dimension*dimension matrix, but need to modify the code in "onestepSGHMC()". 
C         = 3                # tuning parameter, friction. 
learning_rate = 0.001
temperature = 1


## ========= below is the initialization for minibatch sampling narrow potential ravine example ======

# Data = read.csv("dataset/datasetNLR5.csv")
# source("potential/nonlinearReg.R")
# # potential = nonlinearReg
# stochPotential = nonlinearReg
# nData = dim(Data)[1]
# nBatch = max(Data$augment)


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

### =========function 1: SGHMC ===========

SGHMC <- function(temperature, niteration, learning_rate){
  # position:
  x = runif(dimension, min = 0, max = 10)
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


### =========function 2: batchSGHMC===========

batchSGHMC <- function(temperature, nepoch, learning_rate){
  # position:
  # x = runif(dimension, min = 0, max = 1)
  # x = array(x, dim = c(dimension,1))
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # momentum: 
  m = 0
  m = array(m, dim = c(dimension,1))
  
  # result array:
  output = array(NA,dim = c(dimension,2,nepoch*nBatch))   
  output[,1,1] = x      
  output[,2,1] = m
  
  # run SGHMC chain
  for (t in 1:nepoch) {
    for (k in 1: nBatch) {
      Batch = Data[which(Data$augment == k),]
      distribution_sgld <<- function(theta){
        distribution_sgld = stochPotential(theta, Batch, nData)
      }
      temp = onestepSGHMC(x, m, learning_rate, temperature)
      x = temp$position_update
      m = temp$momentum_update
      output[,1,((t-1)*nBatch+k)] = x
      output[,2,((t-1)*nBatch+k)] = m
      
      # print(x)
      # print(t)
    }
  }
  return(output)
}




###=========================================
###=============== run function ============
###=========================================

# output_SGHMC = SGHMC(temperature, 200000, learning_rate)
# 
# plot(output_SGHMC[1,1,1:50000], output_SGHMC[2,1,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGHMC: 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output_SGHMC[1,1,], output_SGHMC[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGHMC: 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")

# write.csv(output_SGHMC, "./variousSGLD/output_SGHMC.csv")


nrun    = 20                  # run nrun independent times.
niteration = 30000*40          # current setup is supposed to >10000
mu_est_SGHMC = array(NA, dim = c(2, nrun))
var_est_SGHMC = array(NA, dim = c(2,2,nrun))

for (i in 1:nrun) {
  print('SGHMC')
  print(paste0(i,"th run"))
  start_time <- Sys.time()
  output_SGHMC = SGHMC(temperature, niteration, learning_rate = 0.001)
  # output_SGHMC = batchSGHMC(temperature, 60000, learning_rate = 0.001)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))
  # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
  # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))


  record = output_SGHMC[,1,c((10000*40):niteration)]

  mu_est_SGHMC[,i] = rowMeans(record)
  var_est_SGHMC[,,i] = cov(t(record))
  print(paste("mu estimate = ", mu_est_SGHMC[,i]))
  print(paste("var estimate = ", var_est_SGHMC[,,i]))
  # print(paste0("1e3 to 2e3: ", rowMeans(record[,1000:2000])))
  # pdf(file=paste0("pRavine",i,".pdf"))
  # plot(record[1,1000:2000],record[2,1000:2000],pch = ".")
  
  # dev.off()
  
  # mu_est_SGHMC[,i] = rowMeans(record[,1000:2000])
  # var_est_EMC[,,i] = cov(t(record))   # mixgaussian
  
  # df = data.frame(mu1 = mu_est_SGHMC[1,], mu2 = mu_est_SGHMC[2,])
  # write.csv(df, "./thetaEstSGHMC.csv")


  
  # print(paste("mu estimate = ", mu_est_SGHMC[,i]))
  # print(paste("var estimate = ", var_est_SGHMC[,,i]))
}



print(paste0('mu est = ',rowMeans(mu_est_SGHMC)))
print(paste0('mu est std = ',sd(mu_est_SGHMC[1,]),', ' ,sd(mu_est_SGHMC[2,])))
print(paste0('var est = ',mean(var_est_SGHMC[1,1,]), ', ' ,mean(var_est_SGHMC[1,2,]), ', ' ,mean(var_est_SGHMC[2,1,]), ', ' ,mean(var_est_SGHMC[2,2,])))
print(paste0('var est std = ',sd(var_est_SGHMC[1,1,]), ', ' ,sd(var_est_SGHMC[1,2,]), ', ' ,sd(var_est_SGHMC[2,1,]), ', ' ,sd(var_est_SGHMC[2,2,])))


##======= plot ========
pdf(file=paste0("mixGaussianSGHMCplot_",i,".pdf"))
for (j in 1:6) {
  # plot(output_snooker$output[1,1,1:(1000*j)], output_snooker$output[2,1,1:(1000*j)], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("EMC ",1000*j, " iteration"), xlab = "x", ylab = "y")
  plot(output_SGHMC[1,1,c(1:(5000*j))*40], output_SGHMC[2,1,c(1:(5000*j))*40], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("SGHMC ",5000*j*40, " iteration"), xlab = "x", ylab = "y")
  points(mu[1,], mu[2,], pch = 3,col = "green")
  
}
# plot(output_snooker$output[1,1,],output_snooker$output[2,1,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "EMC-SGLD sample path", xlab = "x", ylab = "y")
plot(output_SGHMC[1,1,c(10000:30000)*40],output_SGHMC[2,1,c(10000:30000)*40], type = "l", xlim = c(0,10), ylim = c(0,10), main = "SGHMC sample path", xlab = "x", ylab = "y")

dev.off()
# print(paste("avg mu estimate = ", rowMeans(mu_est_SGHMC)))

# mu_est = rowMeans(mu_est_SGHMC)
# sigma_est_1 = mean(var_est_SGHMC[1,1,])
# sigma_est_2 = mean(var_est_SGHMC[2,2,])
# cov_est = mean(var_est_SGHMC[1,2,])
# print(paste0("estimated mu = ", mu_est))
# print(paste0("estimated cov = ", sigma_est_1,", ", sigma_est_2,", ",cov_est))
# 
# df = data.frame(mu1 = mu_est_SGHMC[1,], mu2 = mu_est_SGHMC[2,], var1 = var_est_SGHMC[1,1,], var2 = var_est_SGHMC[2,2,], cov = var_est_SGHMC[1,2,]) # estimating mu and sigma.
# write.csv(df, "./mu_var_est_SGHMC.csv")

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





