### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the SGNHT to produce a chain. algorithm see paper Ding(2014 google research)
### An additional supportive function is defined to deal with the gradient blowup. 

###============= initialization: call functions and initial values=======================
# library(MASS)
# # setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# source("potential/mixtureGaussian.R")
# distribution_sgld <- mixgaussian    

dimension = 2
learning_rate = 0.001
temperature   = 1
A             = 80                    # can be regarded as an initial friction. 

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


### =========function 2: onestepSGNHT===========
## Output one step update through the SGNHT scheme. 

onestepSGNHT <- function(position, momentum, friction, temperature, A, learning_rate = 0.1){
  ## input: 
  # position       : like position of a particle in a 2nd order Langevin Dynamic. dimension*1 array. 
  # momentum       : like momentum of a particle in a 2nd order Langevin Dynamic. dimension*1 array. 
  # friction       : kind of like the friction of the heatbath (not exactly). However it's changing wrt momentum. scaler. 
  # learning_rate  : stepsize of the discretized dynamic.  
  
  dimension = dim(position)[1]
  grad      = distribution_sgld(position)$grad        # calculate the gradient
  friction  = as.numeric(friction)
  wn_length = length(momentum)
  wn        = rnorm(wn_length)
  wn        = array(wn, dim = dim(momentum))            # generate a gaussian noise
  
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
  
  # ================= update ==================
  momentum_update = momentum - friction*momentum*learning_rate - grad*learning_rate + sqrt(2*temperature*A*learning_rate)*wn
  position_update = position + momentum_update*learning_rate
  friction_update = friction + (t(momentum_update)%*%momentum_update/dimension - 1)*learning_rate
  friction_update = min(friction_update, 10)
  
  return(list("position_update" = position_update, "momentum_update" = momentum_update, "friction_update" = friction_update))
}



### ======================================================
### =============== main function=========================
### ======================================================

### =========function 1: SGNHT ===========

SGNHT <- function(A, temperature, niteration, learning_rate){
  # position, dimension*1:
  x = runif(dimension, min = 0, max = 1)
  x = array(x, dim = c(dimension,1))
  # momentum, dimension*1: 
  m = 0
  m = array(m, dim = c(dimension,1))
  # friction, scaler:
  gamma = A
  
  # result array of positions:
  output = array(NA,dim = c(dimension,niteration))   
  output[,1] = x      
  
  # run SGNHT chain
  for (t in 2:niteration) {
    
    temp       = onestepSGNHT(x, m, gamma, temperature, A, learning_rate = 0.1)
    x          = temp$position_update
    m          = temp$momentum_update
    gamma      = temp$friction_update
    output[,t] = x
    
    # print(x)
    # print(t)
  }
  return(output)
}

### =========function 2: batchSGNHT===========

batchSGNHT <- function(A, temperature, nepoch, learning_rate){
  # position, dimension*1:
  # x = runif(dimension, min = 0, max = 1)
  # x = array(x, dim = c(dimension,1))
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # momentum, dimension*1: 
  m = 0
  m = array(m, dim = c(dimension,1))
  # friction, scaler:
  gamma = A
  
  # result array of positions:
  output = array(NA,dim = c(dimension,nepoch*nBatch))  
  output[,1] = x      
  
  # run SGNHT chain
  for (t in 1:nepoch){
    for (k in 1:nBatch){
      Batch = Data[which(Data$augment == k),]
      distribution_sgld <<- function(theta){
        distribution_sgld = stochPotential(theta, Batch, nData)
      }
      temp       = onestepSGNHT(x, m, gamma, temperature, A, learning_rate = 0.1)
      x          = temp$position_update
      m          = temp$momentum_update
      gamma      = temp$friction_update
      output[,t] = x
    }
  }
  
  return(output)
}

###=========================================
###=============== run function ============
###=========================================

# output_SGNHT = SGNHT(A, temperature, 200000, learning_rate)

# plot(output_SGNHT[1,1:50000], output_SGNHT[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGNHT: 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output_SGNHT[1,], output_SGNHT[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGNHT: 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# write.csv(output_SGNHT, "./variousSGLD/output_SGNHT.csv")

nrun    = 10                  # run nrun independent times.
mu_est_SGNHT = array(NA, dim = c(2, nrun))
# var_est_SGNHT = array(NA, dim = c(2,2,nrun))

for (i in 1:nrun) {
  print(paste0(i,"th run"))
  start_time <- Sys.time()
  # output_SGNHT = SGNHT(A, temperature, 150000000, learning_rate = 0.001)
  output_SGNHT = batchSGNHT(A, temperature, 60000, learning_rate = 0.001)
  end_time <- Sys.time()
  print(paste0("running time = ", end_time - start_time))
  # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
  # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))
  
  
  record = output_SGNHT[,c(1:2000)*3000]
  
  # mu_est_SGNHT[,i] = rowMeans(record)
  # var_est_SGNHT[,,i] = cov(t(record))
  print(paste0("1e3 to 2e3: ", rowMeans(record[,1000:2000])))
  pdf(file=paste0("pRavine",i,".pdf"))
  plot(record[1,1000:2000],record[2,1000:2000],pch = ".")
  
  dev.off()
  
  mu_est_SGNHT[,i] = rowMeans(record[,1000:2000])
  # var_est_EMC[,,i] = cov(t(record))   # mixgaussian
  
  df = data.frame(mu1 = mu_est_SGNHT[1,], mu2 = mu_est_SGNHT[2,])
  write.csv(df, "./thetaEstSGNHT.csv")
  
}
#   pdf(file=paste0("mixGaussianSGnhtplot_",i,".pdf"))
#   for (j in 1:50) {
#     plot(record[1,c(1:(1000*j))], record[2,c(1:(1000*j))], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("SGNHT ",1000*j*3000, " iteration"), xlab = "x", ylab = "y")
#     points(mu[1,], mu[2,], pch = 3,col = "green")
#     
#   }
#   plot(record[1,],record[2,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "SGLD sample path", xlab = "x", ylab = "y")
#   
#   dev.off() 
#   print(paste("mu estimate = ", mu_est_SGNHT[,i]))
#   print(paste("var estimate = ", var_est_SGNHT[,,i]))
# }
# 
# mu_est = rowMeans(mu_est_SGNHT)
# sigma_est_1 = mean(var_est_SGNHT[1,1,])
# sigma_est_2 = mean(var_est_SGNHT[2,2,])
# cov_est = mean(var_est_SGNHT[1,2,])
# print(paste0("estimated mu = ", mu_est))
# print(paste0("estimated cov = ", sigma_est_1,", ", sigma_est_2,", ",cov_est))
# 
# df = data.frame(mu1 = mu_est_SGNHT[1,], mu2 = mu_est_SGNHT[2,], var1 = var_est_SGNHT[1,1,], var2 = var_est_SGNHT[2,2,], cov = var_est_SGNHT[1,2,]) # estimating mu and sigma. 
# write.csv(df, "./mu_var_est_SGNHT.csv")

##====================================
## =========== test ==================
##====================================

### test for SGNHT applied on 

# position = mvrnorm(1, mu = rep(0,dimension), Sigma = diag(1,dimension))
# position = array(position,dim = c(2,1))
# momentum = array(c(0,0),dim = c(2,1))
# friction = 1
## x = array(c(1,-0.5,0.3,0.2,0.4),dim = c(1,5))
# learning_rate = 0.1
# onestepSGNHT(position, momentum, friction, A, learning_rate = 0.1)



