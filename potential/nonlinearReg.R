### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation of gradient and density of Nonlinear regression in Kim 2020 section 4.3. 



### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
## parameter setting. 
nSet = 1

nonlinearReg <-function(theta, Batch, nData){
  ## This function calculates the LEVERAGED posterior density and LEVERAGED gradient (-d log pd). 
  ## The leveraged posterior is exp(-U/c) where c is some constant. c is needed because U is usually too large. 
  ## same with the leveraged gradient. 
  ## c here is nData. 
  ## Input: 
  #    theta : dimension of 2*n. 2 = dimension of space, n = number of chains. 
  #    Batch : output[which(output$augment == 2),]
  #    nData : number of data in total dataset.  
  ## Output: 
  #    $pd   : leveraged posterior density
  #    $gradU: leveraged stochastic gradient
  nBatch = dim(Batch)[1]
  nTheta = dim(theta)[1]
  x = matrix(Batch$x)%*%rep(1,nTheta)
  y = matrix(Batch$y)%*%rep(1,nTheta)
  theta1 = t(matrix(theta[,1])%*%rep(1,nBatch))
  theta2 = t(matrix(theta[,2])%*%rep(1,nBatch))
  
  fx = (x - 1)^2 + 2*sin(x*theta1) + theta1/30 + cos(x*theta2 - 1) - theta2/20
  negLogPrior = (theta[,1]^2 + theta[,1]^2)/2
  negStochLogLike = nData/nBatch* colSums((y - fx)^2)/2
  pd = exp(-(negLogPrior+negStochLogLike)/nData)       # nData: leverage
  
  DnegLogPrior = theta
  DnegStochLogLike = matrix(c(colSums(-(y-fx)*(2*x*cos(theta1*x) +1/30)), colSums(theta2 + (y - fx)*(x*sin(theta2*x -1) +1/20)))*nData/nBatch, nrow = dim(theta)[1], ncol = dim(theta)[2])
  gradU = (DnegLogPrior + DnegStochLogLike)/nData               # nData: leverage
  return(list("pd" = pd, "gradU" = gradU))
}


### ====================================================
### ================!!!   TEST    !!!===================
### ====================================================
theta1 = seq(-20,30,0.1)
theta2 = seq(-20,30,0.1)
theta = expand.grid(theta1,theta2)



