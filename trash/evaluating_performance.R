### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the functions for evaluating the different performances of the algorithms, and the block running for result evaluations.


## ========================= ancillary function1 =====================
l2bycol <- function(vec){
  # input:
  #      vec :                 an array of dim = c(40,1), from the columns of diffX_MU (40 * somelength)
  # output:
  #      distances:            a vector of length 20, the l2 distance of a certain position to the 20 mu's

  # first conform vec to dim = c(2,10)
  vec = array(vec,dim = c(2,20))
  vec2 = vec^2
  l2   = sqrt(colSums(vec2))
  return(l2)
}

## ========================= acillary function 2 =============================
## define a function that returns the #of entry that has the min of a vector.

min_index <- function(vec){
  # input:
  #     vec :                   an aarray of dim = c(20,1), from the columns of l2X_MU
  # output:
  #     ind :                   a number between (and contain) 1 to 20, index of the min.
  ind = which(vec == min(vec))
  return(ind)
}

## ========================= acillary function 3 =============================

# lastly define a function that uses a for loop to calculate the sample mean of each class (corresponding to a mu)

calculate_sample_mean <- function(fullindex, mcmcoutput){
  # input:
  # fullindex :                      a vector of length as the output of mcmc, recording for every iteration the class it should belong to.
  # mcmcoutput :                     the output of mcmc, 2*niteration matrix.

  # output:
  # estimated_mu :                   a matrix of dimension 2*20 that are the sample means of the classes.
  # estimated_var:                   a matrix of dimension 2*2, covariance matrix of the within class covariance.
  samplemean = c()
  output_tocenter = mcmcoutput
  for (class in 1:20) {
    classindex = which(fullindex == class)
    if (length(classindex) == 1){
      class_samplemean = mcmcoutput[,classindex]
    }
    else{
      class_samplemean = rowMeans(mcmcoutput[,classindex])
    }

    outputcentered = mcmcoutput[,classindex] - class_samplemean

    output_tocenter[,classindex] = outputcentered
    samplemean = cbind(samplemean,class_samplemean)
  }
  variance_est = cov(t(output_tocenter))

  return(list("estimated_mu" = samplemean, "estimated_var" = variance_est))
}


### ========================================================
### ================!!! MAIN FUNCTION !!!===================
### ========================================================

eval_mu_var <- function(chain, mu){

  # input:
  # chain  :            the markov chain for estimating mu and sigma. 2*somelength.
  # mu     :            underlying mu. 2*20.

  # output:
  # error_mu         : a vector of somelength. evaluates the goodness of the estimation of in-class mu's.
  # estimation_var   : a vector of somelength. the within-group estimation of the variance = sigma^2.

  ## preprocessing: 
  chainlength = dim(chain)[2]

  # conform chain to a matrix of size 40*chainlength by duplicating and attaching by row.
  conform_chain = do.call(rbind, replicate(20, chain, simplify=FALSE))
  # conform mu to a matrix of size 40*chainlength by first conforming to a column vector then duplicating and attaching by column.
  conform_mu = array(mu,dim = c(40,1))
  conform_mu = do.call(cbind, replicate(chainlength, conform_mu, simplify=FALSE))

  ## calculations: 
  # difference by dimension
  diffX_MU = conform_chain - conform_mu              # each row is an iteration.
  # calculate the l2 distance of each sample to the 20 mu's.
  l2X_MU = apply(diffX_MU, 2, l2bycol)              # a 20*chainlength matrix.
  # now calculate which is the closest mu from the positions of the chain.
  min_ind_XtoMU = apply(l2X_MU,2,min_index)

  ## Now evaluate how good the chain estimates mu and variance along the time.
  # the estimation of mu is evaluated by the largest of the 20 in-class l2 distance of the sample means to the true means (at every iteration). It should be smaller the better.
  # the estimation of variance is evaluated by the within-group sample variance (at every iteration). It should be as close as possible to sigma^2 = 0.01.

  # first define two empty vectors to record the evaluations along the iterations.
  error_mu_path = c()
  estimation_var_path = c()

  for (t in seq(1000,chainlength,1000)) {     # start after burn- in period
    # calculate the estimated mean and variance until time = t.
    sample_mu_var = calculate_sample_mean(min_ind_XtoMU[1:t], chain[,1:t])
    sample_mu = sample_mu_var$estimated_mu           # a 2*20 matrix. Might contain NA.
    sample_var = sample_mu_var$estimated_var          # a 2*2 matrix. within group cov matrix of 2 dimensions.

    # define the error in mu to be the largest of the 20 in-class l2 distances of the sample means to the true means.
    error_mu = sqrt(colSums((sample_mu - mu)^2))         # calculate the l2 distances of the sample means to the true means for the 20 classes.
    error_mu = max(error_mu)                         # and take that larges distance. Returns a scaler, could be NA.

    # define the estimation in variance to be the average of the two variances in sample_var.
    estimation_var = mean(diag(sample_var))         # returns a scaler.

    # record the evaluations.
    error_mu_path = c(error_mu_path, error_mu)
    estimation_var_path = c(estimation_var_path, estimation_var)
  }

  return(list("error_mu" = error_mu_path, "estimation_var" = estimation_var_path))

}


### ============================================
### ================run block===================
### ============================================

output_SGLD = read.csv(file = "./variousSGLD/output_SGLD.csv")
output_SGHMC = read.csv(file = "./variousSGLD/output_SGHMC.csv")
output_SGNHT = read.csv(file = "./variousSGLD/output_SGNHT.csv")
output_pSGLD = read.csv(file = "./variousSGLD/output_pSGLD.csv")
output_MSGLD = read.csv(file = "./variousSGLD/output_MSGLD.csv")
output_RESGLD = read.csv(file = "./variousSGLD/output_RESGLD.csv")
output_EMCSGLD = read.csv(file = "./variousSGLD/output_EMCSGLD.csv")

chainlength = 200000

# ===================
chain = as.matrix(output_SGLD[,1:chainlength])
evaluation = eval_mu_var(chain,mu)
eval_SGLD = evaluation

# ===================
output = as.matrix(output_SGHMC[,2:dim(output_RESGLD)[2]])
output = array(output,dim = c(2,2,200000))       # conform the output to dimension*2*niteration. 
chain = output[,1,1:chainlength]
evaluation = eval_mu_var(chain,mu)
eval_SGHMC = evaluation

# ===================
chain = as.matrix(output_SGNHT[,1:chainlength])
evaluation = eval_mu_var(chain,mu)
eval_SGNHT = evaluation

# ===================
chain = as.matrix(output_pSGLD[,1:chainlength])
evaluation = eval_mu_var(chain,mu)
eval_pSGLD = evaluation

# ===================
chain = as.matrix(output_MSGLD[,1:chainlength])
evaluation = eval_mu_var(chain,mu)
eval_MSGLD = evaluation

# ===================
output = as.matrix(output_RESGLD[,2:dim(output_RESGLD)[2]])
output = array(output,dim = c(2,2,200000))       # conform the output to dimension*2*niteration. 
chain = output[,2,1:chainlength]

evaluation = eval_mu_var(chain,mu)
eval_RESGLD = evaluation


# ===================
chain = as.matrix(output_EMCSGLD[,1:chainlength])
evaluation = eval_mu_var(chain,mu)
eval_EMCSGLD = evaluation


# ==== x axis ====================
iteration = seq(1000,chainlength,1000)


## ==========plotting ==========================
# plot(error_mu[7920:length(error_mu)]~iteration, type = "l")
plot(eval_SGHMC$error_mu~iteration, type = "l", col = "red", ylim = c(0,1.5), xlab = "iteration", ylab = "maximum bias in mu")
lines(eval_SGLD$error_mu~iteration, type = "l", col = "dark red")
lines(eval_SGNHT$error_mu~iteration, type = "l", col = "yellow")
lines(eval_pSGLD$error_mu~iteration, type = "l", col = "blue")
lines(eval_MSGLD$error_mu~iteration, type = "l", col = "purple")
lines(eval_RESGLD$error_mu~iteration, type = "l", col = "green")
lines(eval_EMCSGLD$error_mu~iteration, type = "l")
legend("topright", legend=c("SGHMC", "SGLD", "SGNHT", "pSGLD", "MSGLD", "RESGLD","EMCSGLD"), col=c("red", "dark red", "yellow", "blue", "purple", "green", "black"), bty = "n", lwd = 1, horiz = "F",cex = 0.5)


plot(eval_SGHMC$estimation_var~iteration, type = "l", col = "red", ylim = c(0,0.3), xlab = "iteration", ylab = "estimated variance")
lines(eval_SGLD$estimation_var~iteration, type = "l", col = "dark red")
lines(eval_SGNHT$estimation_var~iteration, type = "l", col = "yellow")
lines(eval_pSGLD$estimation_var~iteration, type = "l", col = "blue")
lines(eval_MSGLD$estimation_var~iteration, type = "l", col = "purple")
lines(eval_RESGLD$estimation_var~iteration, type = "l", col = "green")
lines(eval_EMCSGLD$estimation_var~iteration, type = "l")
legend("topright", legend=c("SGHMC", "SGLD", "SGNHT", "pSGLD", "MSGLD", "RESGLD","EMCSGLD"), col=c("red", "dark red", "yellow", "blue", "purple", "green", "black"), bty = "n", lwd = 1, horiz = "F",cex = 0.5)


## ================ record the results ======================
error_mu = cbind(eval_SGHMC$error_mu, eval_SGLD$error_mu, eval_SGNHT$error_mu, eval_pSGLD$error_mu, eval_MSGLD$error_mu, eval_RESGLD$error_mu, eval_EMCSGLD$error_mu)
write.csv(error_mu, "./variousSGLD/error_mu.csv")

estimation_var = cbind(eval_SGHMC$estimation_var, eval_SGLD$estimation_var, eval_SGNHT$estimation_var, eval_pSGLD$estimation_var, eval_MSGLD$estimation_var, eval_RESGLD$estimation_var, eval_EMCSGLD$estimation_var)
write.csv(estimation_var, "./variousSGLD/estimation_var.csv")

