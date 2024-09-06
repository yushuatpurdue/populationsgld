### This file is created by Yushu Huang for population SGMCMC experiments.

## This is the function that runs population SGLD with snooker. 

## input: 
#     Ts                 : temperature ladder. vector. 
#     temperature_select : select temperature in snooker. scaler. 
#     niteration         : number of iterations. scaler. 
#     nchain             : number of the chains in the population to run. scaler. 

## output: 
#     a population of paths. array of dimension = c(dim(state space), nchain, niteration)

###==============================================
###=============== !!!MAIN FUNCTION!!! ==========
###==============================================
library(MASS)
# setwd("/Users/yushuhuang/Documents/pu/research/populationsgld")
### =========== parameters ===================
## ESGLD:
nchain  = 5
# Ts      = 1
Ts      = seq(1,10,9/(nchain-1))
# mutation_rate = 1       # proportion of chains to go through SGLD instead of snooker.
mutation_rate = 0.7        # number of chains to go through SGLD instead of snooker.
# learning_rate = 0.01
learning_rate = seq(0.001,0.021,0.02/(nchain - 1))
# learning_rate = seq(0.01,0.03,0.02/(nchain - 1))

## SGLD: 
# nchain  = 1
# Ts      = rep(1,nchain)
# mutation_rate =  1
# learning_rate = rep(0.001,nchain)

## setting for AB polymer 2d. 
# AB = 'ABBABBABABBAB'
# learning_rate = seq(0.005,0.025,0.02/(nchain - 1))

###============= source function ================
# source('potential/ABpolymer.R')
# potential <- function(x){
#   temp <- function(x){
#     temp = ABpolymer(x,AB)
#   }
#   potential = apply(x, 2, temp)
#   len = length(potential)
#   result = array(unlist(potential),dim = c((length(x)/len +1), len))
#   pd = result[(length(x)/len +1),]
#   grad = result[1:(length(x)/len),]
#   return(list('pd' = exp(-pd),'grad' = grad))
# }
# distribution_sgld <- potential

source("potential/mixtureGaussian.R")
distribution_sgld <- mixgaussian               # target posterior distribution.
potential = mixgaussian

# source("potential/2gaussian.R")
# distribution_sgld <- twogaussian               # target posterior distribution.
# potential = twogaussian

snookerSGLD <- function(Ts, temperature_select, niteration, nchain, learning_rate){
###=========== set path ================
  # if(getwd() != "/Users/yushuhuang/Documents/pu/2020fall/populationsgld"){
  #   setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  # }
  
###============= initialization: call functions and initial values=======================
  ## source the potentials to be used for sampling, then define the distribution accordingly. 
  
  # source("potential/2gaussian.R")
  # distribution_sgld <- twogaussian               # target posterior distribution.  
  
  ## sample initial position x. need to change this according to the target distribution!!
  # temp <- mvrnorm(nchain, mu = rep(2.5,dimension), Sigma = diag(25,dimension))      # initial states
  # x <<- t(temp)                                  # DEFINE AS GLOBAL! or snooker(x, Index) wouldn't be able to read it.
  temp <- runif(nchain*dimension, min = 0, max = 10)
  x <<- matrix(temp, nrow = dimension)
  # print(dim(x))
  source("sampler/snooker.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,niteration))   # result array
  # print((dim(output)))
  output[,,1] = x                               # record initialization
  anchorInd = array(NA,dim = niteration)
  anchor = array(NA,dim = c(dimension,niteration))
  
  
  for (t in 2:niteration) {
    full_index = 1:nchain
    if (mutation_rate>0 && mutation_rate<1){
      Index = sample.int(nchain, ceiling(nchain*(1 - mutation_rate)))        # sample (1 - mutation_rate) portion of the chains to do crossover. 
      # print(Index)
      sgld_index = full_index[-Index]
      
      
      temp = snooker(x, temperature_select, Index, tau = Ts)
      x[,Index] = temp$snooker_update
      
      
      x[,sgld_index] = SGLD(x[,-Index], distribution_sgld, epsilon = learning_rate[sgld_index], tau = Ts[sgld_index])
      
      output[,,t] = x
      anchorInd[t] = temp$anchor1Index
      anchor[,t] = temp$anchor1
    }
    else if (mutation_rate == 1){
      x = SGLD(x, distribution_sgld, epsilon = learning_rate, tau = Ts)
      output[,,t] = x
    }
    
    # print(t)
    # print(x)

  }
  
  if(mutation_rate>0 && mutation_rate<1){
    return(list("output" = output, 'anchorInd' = anchorInd, 'anchor' = anchor))
  }
  else if (mutation_rate == 1){
    return(list("output" = output))
  }
  
  
}


###=========================================
###=============== run function ============
###=========================================

# print(paste0("mutation rate = ", mutation_rate))
# print(paste0("learning rate = ", learning_rate[1], "to ", learning_rate[nchain]))
# 
# 
# ##================ lowest potential of AB polymer ==================
# print("ABpolymer 2d")
# print(paste0("T = ",Ts, ", learning_rate = ", learning_rate))
# print(paste0("mutation rate = ", mutation_rate, ", nchain = ", nchain))
# start_time <- Sys.time()
# output = snookerSGLD(Ts = Ts, temperature_select =0.5, niteration=10000, nchain=nchain, learning_rate = learning_rate)
# end_time <- Sys.time()
# print(paste0("running time = ", end_time - start_time))
# 
# pd_history = potential(output$output[,1,])$pd   # larger the better.
# print(paste0("best pd = ",max(pd_history)))
# write.csv(potential(output$output[,1,])$pd, "./13potential.csv")
# write.csv(output$output[,1,], "./13traj.csv")
# ABpolymer(output$output[,1,200],AB)$pd

##================ run option1: evaluate the mean & var estimation ==================
nrun    = 20                  # run nrun independent times.
niteration = 30000          # current setup is supposed to >10000
mu_est_EMC = array(NA, dim = c(2, nrun))   # mixgaussian
var_est_EMC = array(NA, dim = c(2,2,nrun))   # mixgaussian

for (i in 1:nrun) {
  start_time <- Sys.time()
  output_snooker = snookerSGLD(Ts = Ts, temperature_select = 0.5, niteration, nchain = nchain, learning_rate = learning_rate)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))

  record = output_snooker$output[,1,c((10001):niteration)]   # mixgaussian
  # record = output_snooker$output[,,c((10000*40):niteration)]

  mu_est_EMC[,i] = rowMeans(record)   # mixgaussian
  var_est_EMC[,,i] = cov(t(record))   # mixgaussian
#   
  print(paste("mu estimate = ", mu_est_EMC[,i]))
  print(paste("var estimate = ", var_est_EMC[,,i]))

  ##======= plot ========
  # pdf(file=paste0("2gaussianEMCplot_",i,".pdf"))
#   pdf(file=paste0("2gaussianSGLDplot_",i,".pdf"))
#   for (j in 1:20) {
#     hist(as.numeric(output_snooker$output[1,1,1:(500000*j)]),breaks = seq(-8,13,0.5), freq=FALSE, main = paste0("EMC-SGLD ", 500000*j, " iterations"), xlab = "axis 1")
#   }
#   plot(seq(22/2000,22,22/2000),output_snooker$output[1,1,c(1:800)*1200], type = "l", ylim = c(-2,8), main = "EMC-SGLD sample path", xlab = "CPU time (mins)", ylab = "position")
#
#   dev.off()
}
write.csv(output_snooker$output[,1,], "./output_snooker.csv")





#
print(paste0('mu est = ',rowMeans(mu_est_EMC)))
print(paste0('mu est std = ',sd(mu_est_EMC[1,]),', ' ,sd(mu_est_EMC[2,])))
print(paste0('var est = ',mean(var_est_EMC[1,1,]), ', ' ,mean(var_est_EMC[1,2,]), ', ' ,mean(var_est_EMC[2,1,]), ', ' ,mean(var_est_EMC[2,2,])))
print(paste0('var est std = ',sd(var_est_EMC[1,1,]), ', ' ,sd(var_est_EMC[1,2,]), ', ' ,sd(var_est_EMC[2,1,]), ', ' ,sd(var_est_EMC[2,2,])))
#
# jobid = Sys.getenv("SLURM_JOBID")
# write.csv(mu_est_EMC, paste0("./mu_est_EMC",jobid,".csv"))
# write.csv(var_est_EMC, paste0("./var_est_EMC",jobid,".csv"))
# 

#
##==================== run option2: run chain and plot =================

# # 
pdf(file=paste0("mixGaussianEMCplot_",i,".pdf"))
# pdf(file=paste0("mixGaussianSGLDplot_",i,".pdf"))
    for (j in 1:6) {
      plot(output_snooker$output[1,1,1:(5000*j)], output_snooker$output[2,1,1:(5000*j)], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("EMC ",5000*j, " iteration"), xlab = "x", ylab = "y")
      # plot(output_snooker$output[1,1,c(1:(5000*j))*40], output_snooker$output[2,1,c(1:(5000*j))*40], pch = ".", xlim = c(0,10), ylim = c(0,10), main = paste0("SGLD ",5000*j*40, " iteration"), xlab = "x", ylab = "y")
      points(mu[1,], mu[2,], pch = 3,col = "green")

    }
plot(output_snooker$output[1,1,10000:niteration],output_snooker$output[2,1,10000:niteration], type = "l", xlim = c(0,10), ylim = c(0,10), main = "EMC-SGLD sample path", xlab = "x", ylab = "y")
# plot(output_snooker$output[1,1,c(10000:30000)*40],output_snooker$output[2,1,c(10000:30000)*40], type = "l", xlim = c(0,10), ylim = c(0,10), main = "SGLD sample path", xlab = "x", ylab = "y")
points(mu[1,], mu[2,], pch = 3,col = "green")

dev.off()
  

##==================== mixgaussian: snooker trasition mode =================
# x0 = record[,2303]     # sample before transit. change the index before running code
# x1 = record[,2304]     # sample after transit. change the index before running code
# y = output_snooker$anchor[,12304]     # anchor point. change the index before running code
# 
# plot(x0[1], x0[2], col = "black", pch = 20, xlim = c(0,10), ylim = c(0,10), xlab = "X", ylab = "Y", main = "Mode transition with snooker", type = "n")
# 
# # Add points for x0, x1, and y
# points(mu[1,], mu[2,], pch = 3,col = "green")
# # points(record[1,151:180], record[2,151:180], col = "black", pch = 20)
# # points(x1[1], x1[2], col = "black", pch = 20)
# points(x0[1], x0[2], col = "black", pch = 20)
# points(x1[1], x1[2], col = "black", pch = 20)
# points(y[1], y[2], col = "red", pch = 20)
# 
# # Add arrow from x0 to x1
# arrows(x0[1], x0[2], x1[1], x1[2], col = "blue", length = 0.1, angle = 30, code = 2)

# x0 = array(c(c(1.22926, 1.15905),c(1.3458175, 0.8688436),c(5.929494, 3.922557), c(6.12693, 7.39018), c(6.642096, 6.043990), c(7.241896, 4.052879), c(6.267951, 6.858978), c(6.198575, 6.958819)), dim = c(2,8))
# x1 = array(c(c(3.010722, 2.669501), c(0.4503561, 5.5249136), c(4.144714, 7.304190), c(2.742609, 6.694822), c(7.333233, 3.992747), c(6.970024, 8.186827), c(8.633237, 6.231087), c(1.107694, 8.902122)), dim = c(2,8))
# y = array(c(c(2.709536, 2.414134), c(1.142171, 1.927733), c(4.416181, 6.789840), c(3.000080, 6.747723), c(6.393887, 6.780655), c(6.945630, 8.557753), c(9.152201, 6.093323), c(2.097394, 8.524332)), dim = c(2,8))
# 
# pdf(file=paste0("ModeTransit.pdf"))
# for (i in c(1:8)){
#   plot(x0[1,i], x0[2,i], col = "black", pch = 20, xlim = c(0,10), ylim = c(0,10), xlab = "X", ylab = "Y", main = "Mode transition with snooker", type = "n")
# 
#   # Add points for x0, x1, and y
#   points(mu[1,], mu[2,], pch = 3,col = "green")
#   # points(record[1,151:180], record[2,151:180], col = "black", pch = 20)
#   # points(x1[1], x1[2], col = "black", pch = 20)
#   points(x0[1,i], x0[2,i], col = "black", pch = 20, cex = 3)
#   points(x1[1,i], x1[2,i], col = "black", pch = 20, cex = 3)
#   points(y[1,i], y[2,i], col = "red", pch = 20, cex = 3)
# 
#   # Add arrow from x0 to x1
#   arrows(x0[1,i], x0[2,i], x1[1,i], x1[2,i], col = "blue", length = 0.1, angle = 30, code = 2)
# }
# dev.off()


## ======= mixgaussian below ===========
# mu_est = rowMeans(mu_est_EMC)
# sigma_est_1 = mean(var_est_EMC[1,1,])
# sigma_est_2 = mean(var_est_EMC[2,2,])
# cov_est = mean(var_est_EMC[1,2,])
# print(paste0("estimated mu = ", mu_est))
# print(paste0("estimated cov = ", sigma_est_1,", ", sigma_est_2,", ",cov_est))
# 
# df = data.frame(mu1 = mu_est_EMC[1,], mu2 = mu_est_EMC[2,], var1 = var_est_EMC[1,1,], var2 = var_est_EMC[2,2,], cov = var_est_EMC[1,2,]) # estimating mu and sigma.
# write.csv(df, "./mu_var_est.csv")
## ========= mixgaussian above ============

# # plot(output_snooker$output[1,,], output_snooker$output[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker 20 chains 200000 iteration", xlab = "x", ylab = "y")
# # points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# # 
# plot(output_snooker$output[1,1,], output_snooker$output[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 50000 iteration", xlab = "x", ylab = "y")
#  points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# 
# plot(output_snooker$output[1,1,1:50000], output_snooker$output[2,1,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 50000 iteration", xlab = "x", ylab = "y")
# points(output_snooker$mu[1,], output_snooker$mu[2,], pch = 3,col = "green")
# 
# write.csv(output_snooker$output[1,1,], "./output_snooker.csv")



##===========================
## run just SGLD:
# nchain = 1
# Ts = 0.8
# mutation_rate = 1
# learning_rate = 0.02
# output_SGLD = snookerSGLD(Ts, temperature_select = 0.5, 200000, nchain, learning_rate)
# 
# output = array(output_SGLD$output,dim = c(2,200000))
# output = as.matrix(output)
# 
# plot(output[1,1:50000], output[2,1:50000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGLD 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# plot(output[1,], output[2,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGLD 200000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")
# 
# write.csv(output, "./variousSGLD/output_SGLD.csv")

###=========================================
###=============== unwrapped code ==========
###=========================================


# ###=========== set path ================
# setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
# 
# nchain = 10
# niteration = 5000
# 
# ###============= initialization: call functions and initial values=======================
# source("potential/2gaussian.R")
# distribution_sgld <- twogaussian               # pass potential(r) to a general name that is recognizable for function SGLD.
# x = mvrnorm(nchain, mu = rep(0,dimension), Sigma = diag(1,dimension))
# x = t(x)
# source("sampler/snooker.R")
# 
# ###=============== to be wrapped into a function =====================
# 
# 
# output = array(NA,dim = c(2,10,niteration))
# output[,,1] = x
# 
# for (t in 1:niteration) {
#   full_index = 1:10
#   Index = sample.int(10,5)        # sample 5 integers between 1 and 10 without replacement.
#   sgld_index = full_index[-Index]
# 
#   # distribution_sgld <- snooker_scaler               # pass potential(r) to a general name that is recognizable for function SGLD.
# 
#   x[,Index] = snooker(x,Index)
# 
# 
#   x[,sgld_index] = SGLD(x[,-Index], distribution_sgld)
# 
#   output[,,t] = x
# }







