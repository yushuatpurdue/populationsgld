### This file is created by Yushu Huang for EMC-SGLD experiments.
### It aims to analyze the results returned by the experiments, esp. the final data plotting etc. 

#========================= testing and solving the vanishing gradient =======================
## This block targets at the gradient blow-up in the mixture of gaussian model. 
## setting: 2-dim, sigma = 0.1

## Idea: 
# 1. When a gradient vanishes in a not-too-sparse space, it must be too far away from the centers region, meaning that we need to find the direction of the centers region. 
# 2. When a gradient vanishes in a too-sparse space, it could either be too far away from where the centers are, or it could be amid the centers while too far from any of them, in this case we need to keep moving to any direction until the gradient escapes from being NA. 
# 3. Once moving into the correct direction the sample would be kicked into the centers region, and should be playing around quite well that the vanishing gradient doesn't happen again. 

## Solution: if gradient = NA, then calculate x - mu and sum it up by dimensions, i.e. rowSums(x - mu). If it's positive in a dimension then x needs to be dragged towards negative, so set gradU = 200, else if it's negative then set gradU = -200. This is a bit similar to the reLU idea. 
# library(ggplot2)
# 
# plot(mu[1,], mu[2,])                    # plot the centers
# # x = c(4.3,4.5)
# x = runif(2,0,10)
# x
# calculate_density_grad(mu, sigma, x, weight)
# rowSums(x - mu)
# 
# x = x - sign(rowSums(x - mu))*200*epsilon 

#========================= plotting and analyzing the snookerSGLD outputs =======================
## read data and conform to numeric mode. 
output_snooker4 = read.csv(file = "./output_snooker_4.csv")
output_snooker1 = read.csv(file = "./output_snooker_1.csv")
# output_snooker = read.csv(file = "./output_snooker_3.csv")
# output_snooker_6 = read.csv(file = "./2gaussian/mu_0/lr0p05/output_snooker_6.csv")

output_snooker1 = as.matrix(output_snooker1)
output_snooker1 = array(output_snooker1,dim = c(5,10,2000))       # conform the output to dimension*nchain*niteration. 
# output_snooker_6 = as.matrix(output_snooker_6[,2:dim(output_snooker_6)[2]])
# output_snooker_6 = array(output_snooker_6,dim = c(5,10,5000))

# output_mu = read.csv(file = "./csv_results/0221_mu0p8/mu.csv")
# output_mu = as.matrix(output_mu[,2:dim(output_mu)[2]])

## plot the result
# plot(output_mu[1,], output_mu[2,], pch = ".") 
# hist(as.numeric(output_snooker[1,1:1200]),breaks = seq(-7,12,0.5),freq=FALSE, main = "snooker SGLD")
pdf(file=paste0("mixGaussianSGLDplot_",1,".pdf"))
plot(seq(13/1200,13,13/1200),output_snooker1[1,1,0:1200], type = "l", ylim = c(-2,8), main = "EMC-SGLD sample path", xlab = "CPU time (mins)", ylab = "position")

# scatter plot:  
# plot(output_snooker[1,,], output_snooker[2,,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker 20 chains 5000 iteration", xlab = "x", ylab = "y")
# points(output_mu[1,], output_mu[2,], pch = 3,col = "green")

# ================
# plot(output_snooker[1,1,], output_snooker[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker: chain 1, 10000 iteration", xlab = "x", ylab = "y")
# points(output_mu[1,], output_mu[2,], pch = 3,col = "green")

# plot(output_snooker[1,1,1:3000], output_snooker[2,1,1:3000], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "snooker chain 1, 3000 iteration", xlab = "x", ylab = "y")
# points(output_mu[1,], output_mu[2,], pch = 3,col = "green")

# ================
# ggplot: doesn't seem informative... I prefer the scatter plot. 
# toplot_snooker = data.frame(t(output_snooker[,1,]))
# colnames(toplot_snooker) <- c("x","y")
# ggplot(toplot_snooker, aes(x=x, y=y) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
# 
# toplot_RESGLD = data.frame(t(output_RESGLD$output[,2,]))
# colnames(toplot_RESGLD) <- c("x","y")
# ggplot(toplot_RESGLD, aes(x=x, y=y) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")


#========================= plotting and analyzing other SGLD outputs =======================
## read data and conform to numeric mode. 
# output_SGHMC_test = read.csv(file = "./variousSGLD/output_SGHMC.csv")
# output_SGHMC_test = as.matrix(output_SGHMC_test[,2:dim(output_SGHMC_test)[2]])
# output_SGHMC_test = array(output_SGHMC_test,dim = c(2,2,50000))       # conform the output to dimension*2*niteration. 
# 
# plot(output_SGHMC_test[1,1,], output_SGHMC_test[2,1,], pch = ".", xlim = c(0,10), ylim = c(0,10), main = "SGHMC: 50000 iteration", xlab = "x", ylab = "y")
# points(mu[1,], mu[2,], pch = 3,col = "green")


#==================== evaluating the estimation of mu and variance ==============================

# # conform output_SGNHT to a 40*50000 matrix by duplicating and attaching by row. 
# X = do.call(rbind, replicate(20, output_SGNHT, simplify=FALSE))
# 
# # conform mu to a 40*50000 matrix by first conforming to a column vector then duplicating and attaching by column. 
# MU = array(mu,dim = c(40,1))
# MU = do.call(cbind, replicate(50000, MU, simplify=FALSE))
# 
# # difference by dimension
# diffX_MU = X - MU              # each row is an iteration. 

## define a function that takes a column of diffX_MU and returns 20 l2 distances. 
# l2bycol <- function(vec){
#   # input: 
#   #      vec :                 an array of dim = c(40,1), from the columns of diffX_MU
#   # output: 
#   #      distances:            a vector of length 20, the l2 distance of a certain position to the 20 mu's
#   
#   # first conform vec to dim = c(2,10)
#   vec = array(vec,dim = c(2,20))
#   vec2 = vec^2
#   l2   = sqrt(colSums(vec2))
#   return(l2)
# }

# # calculate the l2 distance of each position to the 20 mu's. 
# l2X_MU = apply(diffX_MU, 2, l2bycol)

## define a function that returns the entry that has the min of a vector. 

# min_index <- function(vec){
#   # input: 
#   #     vec :                   an aarray of dim = c(20,1), from the columns of l2X_MU
#   # output: 
#   #     ind :                   a number between (and contain) 1 to 20, index of the min. 
#   ind = which(vec == min(vec))
#   return(ind)
# }

# # now calculate which is the closest mu from the positions of the chain. 
# min_ind_XtoMU = apply(l2X_MU,2,min_index)

# lastly define a function that uses a for loop to calculate the sample mean of each class (corresponding to a mu)

# calculate_sample_mean <- function(fullindex, mcmcoutput){
#   # input: 
#   # fullindex :                      a vector of length as the output of mcmc, recording for every iteration the class it should belong to. 
#   # mcmcoutput :                     the output of mcmc, 2*niteration matrix. 
#   # output: 
#   # estimated_mu :                   a matrix of dimension 2*20 that are the sample means of the classes.
#   # estimated_var:                   a matrix of dimension 2*2, covariance matrix of the within class covariance. 
#   samplemean = c()
#   output_tocenter = mcmcoutput
#   for (class in 1:20) {
#     classindex = which(fullindex == class)
#     if (length(classindex) == 1){
#       class_samplemean = mcmcoutput[,classindex]
#     }
#     else{
#       class_samplemean = rowMeans(mcmcoutput[,classindex])
#     }
#     
#     outputcentered = mcmcoutput[,classindex] - class_samplemean
#     
#     output_tocenter[,classindex] = outputcentered
#     samplemean = cbind(samplemean,class_samplemean)
#   }
#   variance_est = cov(t(output_tocenter))
#   
#   return(list("estimated_mu" = samplemean, "estimated_var" = variance_est))
# }
# 
# 
# eval_mu_var <- function(chain, mu){
#   
#   # input: 
#   # chain  :            the markov chain for estimating mu and sigma. 2*somelength. 
#   # mu     :            underlying mu. 2*20. 
#   
#   # output: 
#   # error_mu         : a vector of somelength. evaluates the goodness of the estimation of in-class mu's. 
#   # estimation_var   : a vector of somelength. the within-group estimation of the variance = sigma^2. 
#   
#   chainlength = dim(chain)[2]
#   
#   # conform chain to a matrix of size 40*chainlength by duplicating and attaching by row. 
#   conform_chain = do.call(rbind, replicate(20, chain, simplify=FALSE))        
#   # conform mu to a matrix of size 40*chainlength by first conforming to a column vector then duplicating and attaching by column. 
#   conform_mu = array(mu,dim = c(40,1))
#   conform_mu = do.call(cbind, replicate(chainlength, conform_mu, simplify=FALSE))
#   
#   # difference by dimension
#   diffX_MU = conform_chain - conform_mu              # each row is an iteration. 
#   # calculate the l2 distance of each sample to the 20 mu's. 
#   l2X_MU = apply(diffX_MU, 2, l2bycol)              # a 20*chainlength matrix. 
#   # now calculate which is the closest mu from the positions of the chain. 
#   min_ind_XtoMU = apply(l2X_MU,2,min_index)
#   
#   ## Now evaluate how good the chain estimates mu and variance along the time. 
#   # the estimation of mu is evaluated by the largest of the 20 in-class l2 distance of the sample means to the true means (at every iteration). It should be smaller the better. 
#   # the estimation of variance is evaluated by the within-group sample variance (at every iteration). It should be as close as possible to sigma^2 = 0.01. 
#   
#   # first define two empty vectors to record the evaluations along the iterations. 
#   error_mu_path = c()
#   estimation_var_path = c()
#   
#   for (t in seq(1000,chainlength,1000)) {     # start after burn- in period
#     # calculate the estimated mean and variance until time = t. 
#     sample_mu_var = calculate_sample_mean(min_ind_XtoMU[1:t], chain[,1:t])
#     sample_mu = sample_mu_var$estimated_mu           # a 2*20 matrix. Might contain NA.  
#     sample_var = sample_mu_var$estimated_var          # a 2*2 matrix. within group cov matrix of 2 dimensions. 
#     
#     # define the error in mu to be the largest of the 20 in-class l2 distances of the sample means to the true means. 
#     error_mu = sqrt(colSums((sample_mu - mu)^2))         # calculate the l2 distances of the sample means to the true means for the 20 classes. 
#     error_mu = max(error_mu)                         # and take that larges distance. Returns a scaler, could be NA. 
#     
#     # define the estimation in variance to be the average of the two variances in sample_var. 
#     estimation_var = mean(diag(sample_var))         # returns a scaler. 
#     
#     # record the evaluations. 
#     error_mu_path = c(error_mu_path, error_mu)
#     estimation_var_path = c(estimation_var_path, estimation_var)
#   }
#   
#   return(list("error_mu" = error_mu_path, "estimation_var" = estimation_var_path))
#   
# }


## run the evaluation. 
# output_SGLD = read.csv(file = "./variousSGLD/output_SGLD.csv")
# output_SGHMC = read.csv(file = "./variousSGLD/output_SGHMC.csv")
# output_SGNHT = read.csv(file = "./variousSGLD/output_SGNHT.csv")
# output_pSGLD = read.csv(file = "./variousSGLD/output_pSGLD.csv")
# output_MSGLD = read.csv(file = "./variousSGLD/output_MSGLD.csv")
# output_RESGLD = read.csv(file = "./variousSGLD/output_RESGLD.csv")
# output_EMCSGLD = read.csv(file = "./variousSGLD/output_EMCSGLD.csv")
# 
# chain = output_SGNHT[,1:20000]
# chainlength = 20000
# 
# evaluation = eval_mu_var(chain[,1:chainlength],mu)
# eval_SGNHT = evaluation
# 
# error_mu = evaluation$error_mu
# estimation_var = evaluation$estimation_var
# 
# max(which(is.na(evaluation$error_mu)))
# 
# # iteration = (1000 +(7920 - 1)):chainlength
# iteration = seq(1000,chainlength,1000)
# 
# # plot(error_mu[7920:length(error_mu)]~iteration, type = "l")
# plot(error_mu~iteration, type = "l")
# lines(eval_SGNHT$error_mu~iteration, type = "l", col = "green")
# 
# #==================== evaluating the estimation of mu and variance ==============================
# mu_est_EMC = read.csv(file = "./mu_est_EMC_outlier.csv")
# mu_est_EMC = as.matrix(mu_est_EMC)[,2:dim(mu_est_EMC)[2]]
# var_est_EMC = read.csv(file = "./var_est_EMC_outlier.csv")
# var_est_EMC = as.matrix(var_est_EMC)[,2:dim(var_est_EMC)[2]]
# 
# mu_temp = read.csv(file = "./mu_est_EMC10583481.csv")
# mu_temp = as.matrix(mu_temp)[,2:dim(mu_temp)[2]]
# var_temp = read.csv(file = "./var_est_EMC10583481.csv")
# var_temp = as.matrix(var_temp)[,2:dim(var_temp)[2]]
# mu_est_EMC = cbind(mu_est_EMC, mu_temp)
# var_est_EMC = cbind(var_est_EMC,var_temp)
# 
# mu_temp = read.csv(file = "./mu_est_EMC10583483.csv")
# mu_temp = as.matrix(mu_temp)[,2:dim(mu_temp)[2]]
# var_temp = read.csv(file = "./var_est_EMC10583483.csv")
# var_temp = as.matrix(var_temp)[,2:dim(var_temp)[2]]
# mu_est_EMC = cbind(mu_est_EMC, mu_temp)
# var_est_EMC = cbind(var_est_EMC,var_temp)
# 
# mu_temp = read.csv(file = "./mu_est_EMC10611796.csv")
# mu_temp = as.matrix(mu_temp)[,2:dim(mu_temp)[2]]
# var_temp = read.csv(file = "./var_est_EMC10611796.csv")
# var_temp = as.matrix(var_temp)[,2:dim(var_temp)[2]]
# mu_est_EMC = cbind(mu_est_EMC, mu_temp)
# var_est_EMC = cbind(var_est_EMC,var_temp)
# 
# var_est_EMC = array(var_est_EMC, dim = c(2,2,20))



