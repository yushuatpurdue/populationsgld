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
# mu_temp = matrix( c(4.550054, 4.529108, 4.635438, 4.657882, 6.115618, 6.011488, 5.947732, 5.989498), nrow = 2)
# write.csv(mu_est_EMC, "mu_est_EMC_robust.csv")
# 
# var_temp = matrix(c(10.561816, -1.160173, -1.160173,  6.485462, 10.544170, -1.126407, -1.126407,  6.631775, 10.667723, -1.385709, -1.385709,  6.673100, 10.2903327, -0.8232726, -0.8232726,  6.9141576), nrow = 2)
# write.csv(var_est_EMC, "./var_est_EMC_robust.csv")

mu_temp = matrix( c(5.526708, 5.741937, 5.705698, 5.536091, 5.721873, 5.685136, 5.399478, 5.493425, 5.645647, 5.634415, 5.455537, 5.494820,5.509229, 5.418818, 5.547495, 5.570807, 5.696262, 5.605830, 5.701888,  5.604999,    5.748869, 5.768447, 5.831727, 5.565525, 5.727868, 5.781384, 5.655490, 5.859899, 5.655685, 5.734264, 5.501956, 5.710795, 5.803032, 5.805149, 5.599995, 5.833853, 5.825776, 5.726106, 5.863258, 5.519653), nrow = 2)

var_temp = array(c(9.826892, -1.924173, -1.924173, 6.901258, 9.730197, -2.29198, -2.29198,  6.79852, 9.523288, -2.250789, -2.250789,  6.706485, 9.941184, -1.412199, -1.412199,  7.117355, 9.551394, -2.121002, -2.121002, 6.736620, 9.502248, -2.137411, -2.137411, 6.839418, 9.906269, -1.643993, -1.643993, 7.229482, 9.81645, -2.204860, -2.204860, 6.935247, 9.804440, -1.846376, -1.846376, 6.988694, 9.608407, -2.152124, -2.152124, 6.968925, 10.1119970, -0.9103178, -0.9103178, 7.2464439, 9.908575, -1.865083, -1.865083, 7.011583, 9.806334, -2.342979, -2.342979, 6.785476, 9.850799, -1.846480, -1.846480, 6.956723, 9.877427, -1.655223, -1.655223, 7.112698, 9.573518, -2.413355, -2.413355, 6.990674, 9.512486, -2.35991, -2.35991, 6.69419, 9.703246, -1.812056, -1.812056, 6.911862, 9.561753, -2.452352, -2.452352, 6.519380, 9.903049, -1.311377, -1.311377, 6.990044    ), dim =c(2,2,20))

rowMeans(mu_temp)
sd(mu_temp[1,])
sd(mu_temp[2,])
mean(var_temp[1,1,])
mean(var_temp[2,2,])
sd(var_temp[1,1,])
sd(var_temp[2,2,])
mean(var_temp[1,2,])
sd(var_temp[1,2,])


