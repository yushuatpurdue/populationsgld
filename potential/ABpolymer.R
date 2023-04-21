### This file is created by Yushu Huang for EMC-SGLD experiments.
### It contains the function for calculation the Hamiltonian and its gradient of a polymer under AB model
### And several functions to support the calculation. 

#### coding design: 
### Input: 
# X:      Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer. 
# AB:     AB representation of the chain. string of Fibbonacci length. 
# (constraint: N = length(AB). If not, raise an error. )

### Output: 
# H:      Hamiltonian. scaler. 
# gradH:  gradient of H. vector of length d*N. 

### Intermediate steps:
# function 1 to 3 prepares for the calculations of H and gradH. 
# function 4 and 5 computes the two potential terms that adds up to H. 
# function 1: computes the L2distance**2 matrix. 
# function 2: computes (xi+2 - xi+1)T(xi+1 - xi) vector. 
# function 3: computes C(i,j) in V_LJ. 
# function 4: computes V_theta. 
# function 5: computes V_LJ. 


### ================ function 1: distanceSq(x), compute the L2distance**2 matrix ===============
distanceSq <- function(x){
  ## Input: 
  #     x       Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer.
  ## output:
  #     dist2   matrix where dist2ij = (xi - xj)T(xi - xj). Matrix of dimension N*N. 
  
  N     = dim(x)[2]     # N is the number of monomers
  dist2 = array(NA, dim = c(N,N)) 
  ones  = t(rep(1,N))
  ## idea: for every monomer, copy their position as a matrix and take difference to the rest of the monomers only do the upper triangle to save half computation. 
  for (i in 1:N){
    xi = x[,i]%*% ones
    diffx = xi - x[,i:N]
    l2sq = colSums(diffx^2)
    
    dist2[i:N,i] = l2sq
    dist2[i,i:N] = l2sq
    ones = t(ones[-1])
  }
  
  return(dist2)
}

### ================ function 2: neighborProd(x), compute (xi+2 - xi+1)T(xi+1 - xi) vector ===============
neighborProd <- function(x){
  ## Input: 
  #     x       Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer.
  ## output:
  #     result   vector (xi+2 - xi+1)T(xi+1 - xi), length N-2. 
  
  N      = dim(x)[2]     # N is the number of monomers
  temp   = x[,2:N] - x[,1:(N-1)]   # xi+1 - xi
  result = colSums(temp[,2:(N-1)]*temp[,1:(N-2)])
  return(result)
}


### ================ function 3: C(ab), compute the C multiplier in V_LJ ===============
C <-function(ab){
  ## Input: 
  # ab:    a string like "ABABBAB" describing what the polymer chain looks like. 
  ## Output:
  # C:     matrix, size (N-2)*(N-2). lower triangle are all 0. This matrix computes sigma(i,j) with j>=i+2. 
  polymer = as.list(strsplit(ab, "")[[1]])
  N = length(polymer)
  C = array(0, dim = c(N,N)) 
  for(i in 1:N){
    temp1 = polymer[i]
    for(j in i+2:N){
      temp2 = polymer[j]
      if (temp1 == 'A' & temp2 == 'A'){C[i,j] = 1}
      if (temp1 == 'A' & temp2 == 'B'){C[i,j] = -0.5}
      if (temp1 == 'B' & temp2 == 'A'){C[i,j] = -0.5}
      if (temp1 == 'B' & temp2 == 'B'){C[i,j] = 0.5}
    }
  }
  C = C[1:(N-2), 3:N]
  return(C)
}


### ================ function 4: sumVTheta(DX, dx), compute sum(V_theta) in Hamiltonian ===============

sumVTheta <- function(DX, dx){
  ## Input:
  # DX: output of distanceSq(x), the matrix, N*N.
  # dx: output of neighborProd(x), the vector, N-2.
  
  ## Output: 
  # sum(V_theta) term in Hamiltonian. 
  N = dim(DX)[1]
  d1 = diag(DX[1:(N-1), 2:N])      # vector (xi+1 - xi)T(xi+1 - xi). length N-1
  vTheta = (N-2)/4 - sum(dx/sqrt(d1[1:(N-2)]*d1[2:(N-1)]))/4
  return(vTheta)
}

### ================ function 5: sumVLJ(DX, dx), compute sum(V_theta) in Hamiltonian ===============

sumVLJ <- function(DX, Cmatrix){
  ## Input:
  # DX: output of distanceSq(x), the matrix, N*N.
  # Cmatrix: output of C(ab), C multiplier in V_LJ.
  
  ## Output: 
  # sum(V_LJ) term in Hamiltonian. 
  N = dim(DX)[1]
  D2 = DX[1:(N-2), 3:N]
  
  term1 = 1/D2^6
  term1[lower.tri(term1,diag=F)] = 0
  term2 = Cmatrix/D2^3
  result = sum(4*(term1 - term2))
  return(result)
}

### ================ function 6: gradVTheta(X, DX, dx), compute gradient of V_theta  ===============
gradVTheta <- function(x, DX, dx){
  ## Input:
  # x : Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer.
  # DX: output of distanceSq(x), the matrix, N*N.
  # dx: output of neighborProd(x), the vector, N-2.
  
  ## Output: 
  # gradient of V_theta. Matrix of dime (d,N), same dimensionality with dim(X)
  
  ## computation:
  # If you calculate the gradient of V_theta wrt a certain x_i, this x_i could be the x_i of V_theta(i), or the x_{i+1} of V_theta(i-1), or the x_{i+2} of V_theta(i-2). So the gradient is a summation of 7 parts, denoted below. However, i is from 1 to N, same is i-1 and i-2. This makes the gradient of i = 1,2,N-1 and N be different from the rest. To align this difference, I record all the 7 parts as matrices, and add 0 columns to them. 
  
  N = dim(DX)[1]
  d1 = diag(DX[1:(N-1), 2:N])      # vector (xi+1 - xi)T(xi+1 - xi). length N-1
  
  # A B C are the three terms that contributes VTheta. 
  A = 1/sqrt(d1[1:(N-2)])
  C = 1/sqrt(d1[2:(N-1)])
  B = dx
  
  # the following 7 terms adds up to gradient(VTheta). 
  dAidxiBC = array(0, dim = dim(X))
  temp = A^3*(x[,2:(N-1)] - x[,1:(N-2)])
  dAidxiBC[,1:(N-2)] = temp*B*C
  
  dAidxi1BC = array(0, dim = dim(X))
  temp = A^3*(x[,1:(N-2)] - x[,2:(N-1)])
  dAidxi1BC[,2:(N-1)] = temp*B*C
  
  ABdCidxi1 = array(0, dim = dim(X))
  temp = C^3*(x[,3:N] - x[,2:(N-1)])
  ABdCidxi1[,2:(N-1)] = A*B*temp
  
  ABdCidxi2 = array(0, dim = dim(X))
  temp = C^3*(x[,2:(N-1)] - x[,3:N])
  ABdCidxi2[,3:N] = A*B*temp
  
  AdBdxC = array(0, dim = dim(X))
  temp = (x[,2:(N-1)] - x[,3:N])
  AdBdxC[,1:(N-2)] = temp
  
  AdBdx1C = array(0, dim = dim(X))
  temp = x[,1:(N-2)] + x[,3:N] - 2*x[,2:(N-1)]
  AdBdx1C[,2:(N-1)] = temp
  
  AdBdx2C = array(0, dim = dim(X))
  temp = x[,2:(N-1)] - x[,1:(N-2)]
  AdBdx2C[,3:N] = temp
  
  result = -(dAidxiBC + dAidxi1BC + ABdCidxi1 + ABdCidxi2 + AdBdxC + AdBdx1C + AdBdx2C)/4
  return(result)
  
}
### ================ function 7: gradVLJ(X, DX), compute gradient of V_LJ ===============
gradVLJ <- function(x, DX, Cmatrix){
  ## Input:
  # x       : Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer.
  # DX      : output of distanceSq(x), the matrix, N*N.
  # Cmatrix : output of C(ab), C multiplier in V_LJ.
  
  ## Output: 
  # gradient of V_LJ. Matrix of dime (d,N), same dimensionality with dim(X)
  
  ## computation: 
  # This computation makes use of the DX matrix (distance^2 matrix). In order to satisfy that i = 1,...,N-2 and j = i+2,...,N, the super and sub diagonal of DX needs to be set to 0 (derive it and you'll see). 
  
  ## first, modify DX matrix. First divide 1, because DX will appear in the denominator. Then, the diagonal and the super and sub diagonal need to be 0. Name it invDXhollow. 
  d = dim(x)[1]
  N = dim(x)[2]

  invDXhollow = 1/DX
  diag(invDXhollow) <- 0
  invDXhollow[row (invDXhollow) == col (invDXhollow) - 1] <- 0
  invDXhollow[row (invDXhollow) == col (invDXhollow) + 1] <- 0
  
  ## next, setup Chollow from Cmatrix that its diagonal, super-diagonal and subdiagonal are 0. 
  Chollow = array(0, dim = dim(DX))
  Chollow[1:(N-2),3:N] = Cmatrix
  Chollow = Chollow + t(Chollow)
  
  ## calculate gradient column by column (monomer by monomer)
  result = array(NA,dim = dim(x))
  ones  = t(rep(1,N))
  for(i in 1:N){
    temp = x[,i]%*% ones
    result[,i] = rowSums(t( (-24*invDXhollow[i,]^7 + 12*Chollow[i,]*invDXhollow[i,]^4) * t(temp - x) ))
  }
  
  return(result)
}

ABpolymer <- function(X, AB){
  ### Input: 
  # X:      Positions of the monomers in a Fibbonacci AB polymer. matrix of dim (d,N), d space dimension, N length of polymer. 
  # AB:     AB representation of the chain. string of Fibbonacci length. 
  # (constraint: N = length(AB). If not, raise an error. )
  
  ### Output: 
  # H:      Hamiltonian. scaler. 
  # gradH:  gradient of H. vector of length d*N. 
  
  Cmatrix = C(AB)
  DX = distanceSq(X)
  dx = neighborProd(X)
  H = sumVLJ(DX,Cmatrix) + sumVTheta(DX,dx)
  gradH = gradVLJ(X,DX,Cmatrix) + gradVTheta(X,DX,dx)
  result = list("grad" = gradH,"pd" = H)
  return(result)
}


###======== test ==============
AB = 'ABAB'    # access ith monomer by substr(AB, i,i)
Cmatrix = C(AB)

X = t(mvrnorm(4, mu = rep(0,2), Sigma = diag(1,2)))  # d = 2, N = 4, dim(X) = 2,4
DX = distanceSq(X)
dx = neighborProd(X)
sumVLJ(DX,Cmatrix)
sumVTheta(DX,dx)
gradVLJ(X,DX,Cmatrix)
gradVTheta(X,DX,dx)

ABpolymer(X,AB)


