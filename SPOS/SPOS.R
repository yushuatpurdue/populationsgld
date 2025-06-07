### This file is created by Yushu Huang for local SPOS experiments.
### This file runs SPOS with localized K. Algorithm improves Zhang 2020 SPOS paper. 

library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/research/populationsgld/")
source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian

## ================ test =================
# nchain = 10
# dimension = 2
# lr = 0.03 # for now
# Beta = 1  # for now
# eta = 1 # for now
# truncAt = exp(-1)
# niteration = 1000
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# pd_grad = distribution_sgld(x)
# print(pd_grad)

## =============== function: calculate K matrix ====================
# K is an M*M matrix, M = # of particles.
# K = [[k11,k12,k13,...,k1M],[k21,k22,k23,...,k2M],...,[k1M,k2M,k3M,...kMM]]. Note: row and column is flipped against convention. 

kernelComp <- function(x,eta){
  ## Input: 
  #     x       particles.      d*M matrix. d = dimension, M = number of particles.
  #     eta     scaling parameter. scale.
  ## output:
  #     K       kernel matrix. Matrix of dimension M*M. 
  #     gradK   gradient (of second term) of K matrix, summed over the second term. 
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
  gradX = distribution_sgld(x)$grad
  
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
  negStochGrad = - lr*gradX/Beta
  sgldUpdate = x + negStochGrad + sqrt(2*lr/Beta)*wn
  
  ## variational gradient
  K_gradK = kernelComp(x,eta)
  K = K_gradK$K
  gradK = K_gradK$gradK
  VG1 = -gradX%*%t(K)*lr/M 
  VG2 = + gradK*lr/M
  
  xUpdate = sgldUpdate + VG1 + VG2 
  # return(xUpdate)
  return(list("xUpdate" = xUpdate, "negStochGrad" = negStochGrad[,1], "VG1"= VG1[,1], "VG2" = VG2[,1]))
}

## =============== function: run SPOS chains ====================
SPOS <-function(niteration,nchain, learning_rate, eta, Beta){
  x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
  x = t(x)
  output = array(NA,dim = c(dimension,nchain,niteration))   
  output[,,1] = x 
  VG1path = array(NA,dim = c(dimension,niteration-1)) 
  VG2path = array(NA,dim = c(dimension,niteration-1))
  negStochGrad = array(NA,dim = c(dimension,niteration-1))
  for (t in 2:niteration) {
    
    onestep_update  = onestepSPOS(x, learning_rate, eta, Beta)
    x = onestep_update$xUpdate
    output[,,t] = x
    VG1path[,t-1] = onestep_update$VG1
    VG2path[,t-1] = onestep_update$VG2
    negStochGrad[,t-1] = onestep_update$negStochGrad
    
    # print(x)
    # print(t)
  }
  # return(output)
  return(list("output" = output, "negStochGrad" = negStochGrad, "VG1"= VG1path, "VG2" = VG2path))
}


## =============== run experiments =====================

# M =100
niteration = 6000
nrun = 1
Beta = 1
eta = 1
Mlist = c(100)
# Mlist = c(100,300,500,800,1000,1500)

temp = array(NA, dim = c(length(Mlist), 10))

for (i in 1:length(Mlist)){
  M = Mlist[i]
  lr = 0.018
  mu_est_SPOS = array(NA, dim = c(2, nrun))
  var_est_SPOS = array(NA, dim = c(2,2,nrun))
  start_time <- Sys.time()
  for (run in 1:nrun) {
    
    SPOS_result  = SPOS(niteration,M,lr,eta,Beta) # outputSPOS dimension = dim*M*niteration
    outputSPOS  = SPOS_result$output   # outputSPOS dimension = dim*M*niteration
    negStochGrad = SPOS_result$negStochGrad
    VG1SPOS     = SPOS_result$VG1
    VG2SPOS     = SPOS_result$VG2
    
    mu_est_SPOS[,run] = rowMeans(outputSPOS[,1,100:niteration])
    var_est_SPOS[,,run] = cov(t(outputSPOS[,1,100:niteration]))
    
    # pdf(file=paste0("./SPOS/results/4Gaussian/SPOS_M",M,"_beta",Beta,"_",run,".pdf"))
    pdf(file=paste0("./SPOS/results/4Gaussian/SPOS_M",M,"_beta",Beta,"_lr", lr,".pdf"))
    plot(outputSPOS[1,1,c(100:niteration)], outputSPOS[2,1,c(100:niteration)], pch = ".", cex=3, xlim = c(-2,4), ylim = c(-2,4), main = paste0("SPOS ",niteration, " iteration", " M=", M, " lr=", lr, " Beta=", Beta), xlab = "x", ylab = "y")
    points(mu[1,], mu[2,], pch = 3,col = "green")
    
    plot(outputSPOS[1,1,c(100:niteration)], outputSPOS[2,1,c(100:niteration)], type = "n", xlab = "Dimension 1", ylab = "Dimension 2", main = "Trajectory of 2D Samples")
    lines(outputSPOS[1,1,c(100:niteration)], outputSPOS[2,1,c(100:niteration)])
    
    dev.off()
    # pdf(file=paste0("./SPOS/results/4Gaussian/SPOS_M",M,"_beta",Beta,"_",run,"_hist",".pdf"))
    SPOS_hist = hist(outputSPOS[1,1,c(100:niteration)], breaks = seq(-2,4,0.2))
    # dev.off()
    
  }
  end_time <- Sys.time()
  print(paste0("time for sampling ", nrun, "times with M = ", M, " is ", start_time - end_time))
  # df = data.frame(mu1 = mu_est_SPOS[1,], mu2 = mu_est_SPOS[2,], var1 = var_est_SPOS[1,1,], var2 = var_est_SPOS[2,2,], cov = var_est_SPOS[1,2,])
  # write.csv(df, paste0("./thetaEstSPOS_M",M,".csv"))
  temp[i,] = c(rowMeans(mu_est_SPOS), sd(mu_est_SPOS[1,]), sd(mu_est_SPOS[2,]), mean(var_est_SPOS[1,1,]), mean(var_est_SPOS[1,2,]), mean(var_est_SPOS[2,2,]), sd(var_est_SPOS[1,1,]), sd(var_est_SPOS[1,2,]), sd(var_est_SPOS[2,2,]))
  print(temp[i,])
  
}

df = data.frame(temp)
colnames(df) <- c("mu1", "mu2", "sd_mu1", "sd_mu2", "var1", "var2", "cov", "sd_var1", "sd_var2", "sd_cov")
write.csv(df, paste0("./SPOS_lr",lr,"_beta",Beta, ".csv"))


# print(paste0('mu est = ',rowMeans(mu_est_SPOS)))
# print(paste0('mu est std = ',sd(mu_est_SPOS[1,]),', ' ,sd(mu_est_SPOS[2,])))
# print(paste0('var est = ',mean(var_est_SPOS[1,1,]), ', ' ,mean(var_est_SPOS[1,2,]), ', ' ,mean(var_est_SPOS[2,1,]), ', ' ,mean(var_est_SPOS[2,2,])))
# print(paste0('var est std = ',sd(var_est_SPOS[1,1,]), ', ' ,sd(var_est_SPOS[1,2,]), ', ' ,sd(var_est_SPOS[2,1,]), ', ' ,sd(var_est_SPOS[2,2,])))


# print(df)
# 
# print(sum(outputSPOS[,1,]^2)/2000 - 77.08546)

# accOfMSE <-function(M,niteration,nrun){
#   print(paste0('SPOS, M = ',M,', niteration = ', niteration, ', nrun = ', nrun))
#   ## meanMSE calculates the mean estimation of E(X^2) every 500 iterations
#   ## Input: 
#   #   M         : number of particles
#   #   niteration: number of iterations per run
#   #   nrun      : number of runs. 
#   
#   ## Output:
#   #   meanMSE   : vector of length niteration/500. 
#   
#   acc = rep(0, niteration/500)
#   for (run in 1:nrun){
#     outputSPOS = SPOS(niteration,M,0.03,1,1)     # outputSPOS dimension = dim*M*niteration
#     
#     ## If using only 1st chain:
#     reform = array(outputSPOS[,1,],dim = c(2,500,niteration/500))      # take chain 1 and reform it into groups of 500 data. dim(reform) = 2, 500, #groups
#     SS = apply(reform^2,3,sum)       # for each group, calculate sum of square
#     MS = cumsum(SS)/seq(500,niteration,500)      # calculate the mean square of cumulated sum. This gives the estimator of E(X^2) for #iteration = 500, 1000, 1500, ...
#     
#     ## If want to use all M chains:
#     # reform = array(outputSPOS,dim = c(2,500*M,niteration/500))
#     # SS = apply(reform^2,3,sum)       # for each group, calculate sum of square
#     # MS = cumsum(SS)/(seq(500,niteration,500)*M)      # calculate the mean square of cumulated sum. This gives the estimator of E(X^2) for #iteration = 500, 1000, 1500, ...
#     # 
#     temp = abs(MS - 77.08546)        # 77.08546 = mu1^2 + mu2^2 + var1 + var2
#     # 
#     acc = acc + temp
#     # print(paste0('run ',nrun,' MSE = ',temp))
#   }
#   accOfMSE = acc/nrun
#   print(paste0('at iterations ',seq(500,niteration,500),', accuracy = ',accOfMSE))
#   return(accOfMSE)
# }

## ==========run============
# ACC = c()
# for (M in c(100,300,500,800,1000,1500)){
# # for (M in c(1800,2000,2300)){
#   acc = accOfMSE(M,niteration,nrun)
#   ACC = rbind(ACC,acc)
# }
# 
# df = data.frame(M = c(100,300,500,800,1000,1500), ACC = ACC)
# # df = data.frame(M = c(1800,2000,2300), ACC = ACC)
# write.csv(df, "./MSE_SPOS.csv")
