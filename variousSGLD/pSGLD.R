### This file is created by Yushu Huang for EMC-SGLD experiments.
### This file runs the pSGLD to produce a chain. algorithm see paper Ding(2014 google research)


###============= initialization: call functions and initial values=======================
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
source("potential/mixGauss2order.R")
distribution_sgld <- mixgaussian2
# set.seed(1)
dimension     = 2
learning_rate = 0.01
alpha         = 0.01                 # tuning parameter 1, value in [0,1]
lambda        = 0.08                 # tuning parameter 2

## ========= below is the initialization for minibatch sampling narrow potential ravine example ======

# Data = read.csv("dataset/datasetNLR5.csv")
# source("potential/nonlinearReg2order.R")
# potential = nonlinearReg
# stochPotential = nonlinearReg2ord
# nData = dim(Data)[1]
# nBatch = max(Data$augment)


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

### =========function 1: pSGLD ===========

pSGLD <-function(lambda, alpha, niteration, learning_rate){
  
  # position, dimension*1:
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # x = runif(dimension, min = 0, max = 10)
  # x = array(x, dim = c(dimension,1))
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

### =========function 2: batchSGNHT===========
batchpSGLD <-function(lambda, alpha, nepoch, learning_rate){
  
  # position, dimension*1:
  x = rnorm(dimension, 0,1 )
  x = array(x, dim = c(dimension,1))
  # V, dimension*1: 
  V = 0
  V = array(V, dim = c(dimension,1))
  
  # result array of positions:
  output = array(NA,dim = c(dimension,nepoch*nBatch))   
  output[,1] = x      
  
  # run SGHMC chain
  for (t in 1:nepoch){
    for (k in 1:nBatch){
      Batch = Data[which(Data$augment == k),]
      distribution_sgld <<- function(theta){
        distribution_sgld = stochPotential(theta, Batch, nData)
      }
      temp = onesteppSGLD(x, V, alpha, lambda, learning_rate)
      x = temp$position_update
      V = temp$V_update
      
      output[,((t-1)*nBatch+k)] =x
    }
  }
  return(output)
}

###=========================================
###=============== run function ============
###=========================================

# output_pSGLD = pSGLD(lambda, alpha, 200000, learning_rate)
# 
# plot(output_pSGLD[1,1:50000], output_pSGLD[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "pSGLD: 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output_pSGLD[1,], output_pSGLD[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "pSGLD: 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")

# write.csv(output_pSGLD, "./variousSGLD/output_pSGLD.csv")


nrun    = 20                  # run nrun independent times.
mu_est_pSGLD = array(NA, dim = c(2, nrun))
var_est_pSGLD = array(NA, dim = c(2,2,nrun))
niteration = 30000*40

for (i in 1:nrun) {
  print('pSGLD')
  print(paste0(i,"th run"))
  start_time <- Sys.time()
  output_pSGLD = pSGLD(lambda, alpha, niteration = niteration, learning_rate = 0.001)
  # output_pSGLD = batchpSGLD(lambda, alpha, 50000, learning_rate = 0.001)
  end_time <- Sys.time()
  print(paste0("running time = ", end_time - start_time))
  # hist(as.numeric(output_snooker$output[1,1,]),breaks = seq(-8,13,0.5),main = paste0("EMC-SGLD at ", i, "th run"), xlab = "axis 1")
  # write.csv(output_snooker$output, paste0("./output_snooker_",i,".csv"))
  
  
  record = output_pSGLD[,c((10000*40):niteration)]
  mu_est_pSGLD[,i] = rowMeans(record)
  var_est_pSGLD[,,i] = cov(t(record))
  # record = output_pSGLD[,c(1:2000)*2500]
  # print(paste0("1e3 to 2e3: ", rowMeans(record[,1000:2000])))
  # pdf(file=paste0("pRavinepSGLD",i,".pdf"))
  # plot(record[1,1000:2000],record[2,1000:2000],pch = ".")
  # 
  # dev.off()
  
  # mu_est_pSGLD[,i] = rowMeans(record[,1000:2000])
  # var_est_pSGLD[,,i] = cov(t(record))
  
  # df = data.frame(mu1 = mu_est_pSGLD[1,], mu2 = mu_est_pSGLD[2,])
  # write.csv(df, "./thetaEstpSGLD.csv")
  
  # pdf(file=paste0("mixGaussianpSGLDplot_",i,".pdf"))
  # for (j in 1:50) {
  #   plot(record[1,c(1:(1000*j))], record[2,c(1:(1000*j))], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("pSGLD ",1000*j*3000, " iteration"), xlab = "x", ylab = "y")
  #   points(mu[1,], mu[2,], pch = 3,col = "green")
  #   
  # }
  # plot(record[1,],record[2,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "pSGLD sample path", xlab = "x", ylab = "y")
  
  # dev.off() 
  print(paste("mu estimate = ", mu_est_pSGLD[,i]))
  print(paste("var estimate = ", var_est_pSGLD[,,i]))
}

print(paste0('mu est = ',rowMeans(mu_est_pSGLD)))
print(paste0('mu est std = ',sd(mu_est_pSGLD[1,]),', ' ,sd(mu_est_pSGLD[2,])))
print(paste0('var est = ',mean(var_est_pSGLD[1,1,]), ', ' ,mean(var_est_pSGLD[1,2,]), ', ' ,mean(var_est_pSGLD[2,1,]), ', ' ,mean(var_est_pSGLD[2,2,])))
print(paste0('var est std = ',sd(var_est_pSGLD[1,1,]), ', ' ,sd(var_est_pSGLD[1,2,]), ', ' ,sd(var_est_pSGLD[2,1,]), ', ' ,sd(var_est_pSGLD[2,2,])))


pdf(file=paste0("mixGaussianpSGLDplot_",i,".pdf"))
for (j in 1:6) {
  # plot(output_snooker$output[1,1,1:(1000*j)], output_snooker$output[2,1,1:(1000*j)], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("EMC ",1000*j, " iteration"), xlab = "x", ylab = "y")
  plot(output_pSGLD[1,c(1:(5000*j))*40], output_pSGLD[2,c(1:(5000*j))*40], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("pSGLD ",5000*j*40, " iteration"), xlab = "x", ylab = "y")
  points(mu[1,], mu[2,], pch = 3,col = "green")
  
}
# plot(output_snooker$output[1,1,],output_snooker$output[2,1,], type = "l", xlim = c(0,10), ylim = c(0,10), main = "EMC-SGLD sample path", xlab = "x", ylab = "y")
plot(output_pSGLD[1,c(10000:30000)*40],output_pSGLD[2,c(10000:30000)*40], type = "l", xlim = c(0,10), ylim = c(0,10), main = "pSGLD sample path", xlab = "x", ylab = "y")

dev.off()

# print(paste("avg mu estimate = ", rowMeans(mu_est_pSGLD)))

# mu_est = rowMeans(mu_est_pSGLD)
# sigma_est_1 = mean(var_est_pSGLD[1,1,])
# sigma_est_2 = mean(var_est_pSGLD[2,2,])
# cov_est = mean(var_est_pSGLD[1,2,])
# print(paste0("estimated mu = ", mu_est))
# print(paste0("estimated cov = ", sigma_est_1,", ", sigma_est_2,", ",cov_est))

# df = data.frame(mu1 = mu_est_pSGLD[1,], mu2 = mu_est_pSGLD[2,], var1 = var_est_pSGLD[1,1,], var2 = var_est_pSGLD[2,2,], cov = var_est_pSGLD[1,2,]) # estimating mu and sigma. 
# write.csv(df, "./mu_var_est_pSGLD.csv")


##====================================
## =========== test ==================
##====================================

### test for onesteppSGLD applied on 
# V = array(0, dim = c(dimension,1))




