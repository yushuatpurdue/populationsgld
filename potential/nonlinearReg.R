### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation of gradient and density of Nonlinear regression in Kim 2020 4.3. 
### And one additional function to support the calculation. 
### It also contains a function for generating the regression dataset. The 5 datasets only need to be generated once. 

###============= generate dataset ==========

genDataset <- function(nData){
  x = runif(nData,-2,4)   # generate nData number of sample from unif[-2,4] as x
  fx = (x-1)^2 + 2*sin(20*x) +cos(10*x-1) +1/6
  noise = rnorm(nData,0,1)
  y = x+ noise
  return(list('x' = x,'y' = y))
}

