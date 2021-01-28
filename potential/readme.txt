Any file in this folder contains only functions for calculation of certain potential. These functions should have the following design: 


1. main function. It should be named after the distribution. 

input:  
      positions = matrix of p*n. (p = dimension of space, n = number of chains. )
      parameters if cannot avoid. 
       
output: 
      a list of two arguments 
          $grad = gradient of the potential at x. 
          $pd   = (posterior) probability density at x. 

2. parameters. They should be defined within one block, usually right before the main function block. 

3. supportive functions. They should be defined at the beginning of the file so the parameter and main function definition may call them. 



