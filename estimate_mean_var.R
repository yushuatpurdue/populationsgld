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




output_SGLD = snookerSGLD(0.8, temperature_select = 0.5, 200000, nchain, 0.02)
# hist(as.numeric(output_snooker$output[2,,1:10]),breaks = seq(-7,12,0.5),main = "snooker SGLD")
output = output_SGLD

nrun    = 20                  # run nrun independent times. 
mu_est = array(NA, dim = c(2, nrun))
var_est = array(NA, dim = c(2,2,nrun))
for (i in 1:20) {
  output_MSGLD = MSGLD(smoothingfactor, biasfactor, temperature, 200000, learning_rate)
  output = output_MSGLD
  mu_est[,i] = rowMeans(output)
  var_est[,,i] = cov(t(output))
}
rowMeans(mu_est)
sd(mu_est[1,])
sd(mu_est[2,])

mean(var_est[1,1,])
sd(var_est[1,1,])
mean(var_est[2,2,])
sd(var_est[2,2,])
mean(var_est[1,2,])
sd(var_est[1,2,])

##============== test block =====================
mu_temp = matrix( c(4.550054, 4.529108, 4.635438, 4.657882, 6.115618, 6.011488, 5.947732, 5.989498), nrow = 2)
write.csv(mu_est_EMC, "mu_est_EMC_robust.csv")

var_temp = matrix(c(10.561816, -1.160173, -1.160173,  6.485462, 10.544170, -1.126407, -1.126407,  6.631775, 10.667723, -1.385709, -1.385709,  6.673100, 10.2903327, -0.8232726, -0.8232726,  6.9141576), nrow = 2)
write.csv(var_est_EMC, "./var_est_EMC_robust.csv")



