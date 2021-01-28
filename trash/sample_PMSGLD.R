# ===================================
## population momentum sgld sampler


sample_PMSGLD <- function(iteration, x, m)
{
  ## initializing the parameters. 
  rate = 0.1   # learning rate
  beta = 0.9   # decaying rate of m. 
  a = 0.1      # multiplier on m when updating position.
  temperature = 10
  
  ## initializing the path record.  
  X = array(x, dim = c(dimension,nchain,1))
  ## initializing the max density record. 
  PD = c()
  
  for(k in 1:iteration)
  {
    ## first calculate the density and gradient
    calc_density_grad = apply(x, 2, twogaussian_density_grad)
    grad = c()
    pd   = c()
    for(i in c(1:dim(x)[2]))
    {
      grad = cbind(grad, calc_density_grad[i][[1]]$gradU)  # gradient matrix. each column a full gradient for each chain. 
      pd   = cbind(pd, calc_density_grad[i][[1]]$pd)     # vector of probability density
    }
    PD = cbind(PD, pd)           # record the density
    
    ## then find optimal position
    opt_position = x[,which.max(pd)]
    
    
    ## then update m and x. 
    m = beta*m - (1-beta)*(opt_position%*%t(rep(1, dim(x)[2])) - x)  # update m.
    epsilon = rnorm(nchain*dimension, mean = 0, sd = 0.7)
    epsilon = matrix(epsilon, nrow = dimension)
    x = x - (grad + a*m)*rate + sqrt(rate*temperature)*epsilon 
    #    x = x - (grad)*rate + epsilon 
    
    ## then record the position into the path
    X = array(c(X,x),dim = c(dimension,nchain,dim(X)[3]+1))
  }
  
  return(list("X" = X,"PD" = PD))
}



#======================
## test

# sample_path_p = sample_PMSGLD(50000, x,m)
# plot(sample_path_p$X[1,3,],type = "l")
# hist(sample_path_p$X[1,,1:50000],breaks = seq(-11,16,0.5))
# hist_population = hist(sample_path_p$X[1,,],breaks = seq(-11,16,0.5))
# 
# sample_path = sample_PMSGLD(200000, x,m)
# plot(sample_path$X[1,3,],type = "l")
# hist(sample_path$X[1,,1:200000],breaks = seq(-6,11,0.5))
# hist_parallel = hist(sample_path$X[1,,],breaks = seq(-6,11,0.5))
# 
# sample_path_p$X = c(sample_path_p$X, sample_path_p1$X)
# sample_path_p$X = array(sample_path_p$X,dim = c(2,10,200020))
