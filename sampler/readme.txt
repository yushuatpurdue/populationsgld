Any file in this folder contains only functions for one-step sampling of a Markov chain. These functions should have the following design: 


1. main function. It should be named after the sampling scheme. 

input:  
      positions = matrix of p*n. (p = dimension of space, n = number of chains. )
      parameters if cannot avoid. 
       
output: 
      One-step update of the positions inputed. Usually named with format *_update. 

2. Initialization. This should be done within one block, usually at the very beginning of the file. The initialization block may define global parameters and source other relevant files. 

3. supportive functions. They should be defined before the main function so that the main can call them. 



