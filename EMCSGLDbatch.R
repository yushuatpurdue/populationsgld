library(MASS)

Data = read.csv("dataset/datasetNLR1.csv")
source("potential/nonlinearReg.R")
# potential = nonlinearReg
stochPotential = nonlinearReg
nData = dim(Data)[1]
nBatch = max(Data$augment)
## note: need to specify batchs and nData. 

snookerSGLDbatch <- function(Ts, temperature_select, nepoch, nchain, learning_rate){
  ###=========== set path ================
  # if(getwd() != "/Users/yushuhuang/Documents/pu/2020fall/populationsgld"){
  #   setwd("/Users/yushuhuang/Documents/pu/2020fall/populationsgld")
  # }
  
  ###============= initialization: call functions and initial values=======================
  
  theta <<- t(mvrnorm(nchain, mu = c(0,0), Sigma = diag(1,2)))
  
  source("sampler/snookerBatch.R")                    # source one-step snooker sampler. Note: SGLD sampler has got sourced in "sampler/snooker.R", no need to source again. 
  
  ###=============== run population SGLD with snooker =====================
  
  output = array(NA,dim = c(dimension,nchain,nepoch*nBatch))   # result array
  # output[,,1] = x                               # record initialization
  
  for (t in 1:nepoch) {
    for (k in 1:nBatch) {
      Batch = Data[which(Data$augment == k),]
      potential <<- function(theta){
        potential = stochPotential(theta, Batch, nData)
      }
      full_index = 1:nchain
      if (mutation_rate>0 && mutation_rate<1){
        Index = sample.int(nchain, ceiling(nchain*(1 - mutation_rate)))        # sample (1 - mutation_rate) portion of the chains to do crossover. 
        # print(Index)
        sgld_index = full_index[-Index]
        
        
        theta[,Index] = snooker_batch(theta, temperature_select, Index, tau = Ts)
        
        
        theta[,sgld_index] = SGLD(theta[,-Index], distribution_sgld = potential, epsilon = learning_rate[sgld_index], tau = Ts[sgld_index])
      }
      else if (mutation_rate == 1){
        theta = SGLD(theta, distribution_sgld = potential, epsilon = learning_rate, tau = Ts)
      }
      
      
      output[,,((t-1)*nBatch+k)] <- theta
      # print(theta)
    }
    
    
    
    # print(t)
  }
  
  return(list("output" = output))
  
}

## test
# Batch = Data[which(Data$augment == 1),]
# nData = 10000
# potential <- function(theta){
#   temp = nonlinearReg(theta, Batch, nData)
# }
# theta = mvrnorm(20, mu = c(0,0), Sigma = diag(1,2))
# nchain  = 20
# nepoch = 10

## ================ run simulation ======================
nchain  = 20
Ts      = seq(1,5,4/(nchain-1))
# Ts      = rep(1,nchain)
mutation_rate = 0.8        # number of chains to go through SGLD instead of snooker.
learning_rate = seq(0.001,0.021,0.02/(nchain - 1))

# output_snookerBatch = snookerSGLDbatch(Ts, temperature_select = 0.5, nepoch = 40, nchain = 10, learning_rate)

nrun    = 20                  # run nrun independent times.
mu_est_EMC = array(NA, dim = c(2, nrun))   # mixgaussian
# var_est_EMC = array(NA, dim = c(2,2,nrun))   # mixgaussian

for (i in 1:nrun) {
  start_time <- Sys.time()
  output_snookerBatch = snookerSGLDbatch(Ts, temperature_select = 0.5, nepoch = 20, nchain = 20, learning_rate)
  end_time <- Sys.time()
  print(paste0("running time for ", i, "th run = ", end_time - start_time))

  # record = output_snookerBatch$output[,1,]   # mixgaussian
  # print(paste0("1e3 to 2e3: ", rowMeans(output_snookerBatch$output[,1,1000:2000]), "; 2e3 to 3e3: ", rowMeans(output_snookerBatch$output[,1,2000:3000]), "; 3e3 to 4e3: ", rowMeans(output_snookerBatch$output[,1,3000:4000])))
  print(paste0(rowMeans(output_snookerBatch$output[,1,1000:2000])))
  # pdf(file=paste0("pRavine",i,".pdf"))
  # plot(output_snookerBatch$output[1,1,1000:2000],output_snookerBatch$output[2,1,1000:2000],pch = ".")
  # 
  # dev.off()

  mu_est_EMC[,i] = rowMeans(output_snookerBatch$output[,1,1000:2000])
  # var_est_EMC[,,i] = cov(t(record))   # mixgaussian

  df = data.frame(mu1 = mu_est_EMC[1,], mu2 = mu_est_EMC[2,])
  write.csv(df, "./thetaEst.csv")
  
}
