### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation of gradient and density of Nonlinear regression in Kim 2020 4.3. 
### And one additional function to support the calculation. 

###============= calculate_density_grad: supportive function ==========

## calculate_density_grad: a function calculate the posterior density and gradient of - log density wrt a given x
## input:
#    dataSet: (x,y) dataset for nonlinear regression. generated from the genDataset function. a dataframe of 2 arguments: $x and $y
#    x:       The current position of a sample.  

## output:
#    $pd: The posterior density at x. 
#    $gradU: The gradient of -log posterior. 

#======================
calculate_density_grad <- function(dataSet, theta)
{
  
}

### ====================================================
### ================!!! MAIN BODY !!!===================
### ====================================================
## parameter setting. 
nSet = 1
# ===================================
## calculate the probability density and gradient for a position x.
## input: 
#    x:   positions. matrix of p*n. p = dimension of space, n = number of chains. 
## output: 
#    a data frame of two columns
#         $grad = gradient of the potential at x. 
#         $pd   = (posterior) probability density at x. 

## side note: potential = -log(probability density)
nonlinearReg <-function(x){
  
}
