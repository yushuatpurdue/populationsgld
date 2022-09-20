### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains a function for generating the regression dataset. We want 5 datasets that only need to be generated once. 

###============= function for generate dataset ==========
## genDataset is a function that generates nData number of (x,y) pairs from the nonlinear regression:
## y = (x-1)^2 + 2*sin(theta1*x) + theta1/30 + cos(theta2*x-1) - theta2/20 + e ~ N(0,1)
# Input:
#   nData : number of (x,y) datapoints to be generated
#   nSet  : id of the dataset. 
# output:
#   $x    : list of length nData. explanatory variable
#   $y    : list of length nData. response variable

genDataset <- function(nData, nSet, batchSize){
  # generate x and y
  x      = runif(nData,-2,4)                 # generate nData number of sample from unif[-2,4] as x
  fx     = (x-1)^2 + 2*sin(20*x) +cos(10*x-1) +1/6  # underlying true model f(x)
  noise  = rnorm(nData,0,1)
  y      = fx+ noise                 # y = f(x) + noise
  
  # augment dataset. 
  one = rep(1, batchSize)
  nBatches = as.integer(nData/batchSize)+1
  augment = one
  for (i in 2:nBatches){
    temp = one*i
    augment = append(augment,temp)
  }
  augment = augment[1:nData]
  
  output = data.frame('augment' = augment, 'x' = x,'y' = y) # use which(output$augment == 1) to find batch. 
  write.csv(output, paste0("./dataset/datasetNLR",nSet,".csv"), row.names=FALSE)
  
}


###============= run: generating dataset ==========
## nSet datasets are generated, each containing nData datapoints. 

## parameter setting
nSet = 5
nData = 10000
batchSize = 100

for (i in 1:nSet)
{
  if (file.exists(paste0("./dataset/datasetNLR",i,".csv"))){
    read.csv(paste0("./dataset/datasetNLR",i,".csv"))
  } else {
    genDataset(nData,i, batchSize)
  }
}

# ========== test ============
# genDataset(10000,1)
# xy = read.csv("./dataset/datasetNLR1.csv")
# file.exists("./dataset/datasetNLR2.csv")

# x = seq(-2,4,0.01)
# fx = (x-1)^2 + 2*sin(20*x) +cos(10*x-1) +1/6 + rnorm(length(x),0,1)
# plot(x,fx,'l')