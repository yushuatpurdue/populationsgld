### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for estimating the overall sample mean and variance. 

##============== calculate the underlying mean and variance =====================

mu_overall = rowMeans(mu)
var_betweengroup = (mu - mu_overall)%*%t(mu - mu_overall)/20
var_overall = var_betweengroup + diag(0.01,2)

output = as.matrix(output_RESGLD[,2:dim(output_RESGLD)[2]])
output = array(output,dim = c(2,2,200000))       # conform the output to dimension*2*niteration. 
chain = output[,2,1:200000]

output = as.matrix(output_SGHMC[,2:dim(output_RESGLD)[2]])
output = array(output,dim = c(2,2,200000))       # conform the output to dimension*2*niteration. 
chain = output[,1,1:200000]

chain = output_snooker
rowMeans(chain)
cov(t(chain))


output_snooker = read.csv(file = "./variousSGLD/output_snooker.csv")



