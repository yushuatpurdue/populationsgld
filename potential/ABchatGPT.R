# Function to calculate the Lennard-Jones potential between two monomers
#
# @param r: Numeric scalar, the distance between two monomers
# @param sigma: Numeric scalar, the distance at which the potential is zero
# @param c: Numeric scalar, the scaling factor for the potential
#
# @return: Numeric scalar, the Lennard-Jones potential between two monomers
calculate_LJ_potential <- function(r, sigma, c) {
  4 * ((sigma / r)^12 - c * (sigma / r)^6)
}

# Function to calculate the bond angle potential between two consecutive bonds in the polymer
#
# @param bond_angle: Numeric scalar, the angle between two consecutive bonds in the polymer
#
# @return: Numeric scalar, the bond angle potential between two consecutive bonds in the polymer
calculate_bond_potential <- function(bond_angle) {
  K_bond <- 1/4
  K_bond * (1 - cos(bond_angle))
}

# Function to calculate the Hamiltonian and its gradient of a polymer chain
#
# @param positions: Numeric matrix of dimension N x d, where N is the number of monomers in the polymer and d is the dimension of the state space. This matrix records the positions of the N monomers in the polymer.
# @param types: Character vector of length N, recording the types of individual monomers.
#
# @return: List of two elements:
#   1. Numeric scalar, the Hamiltonian of the polymer chain
#   2. Numeric matrix of dimension N x d, the gradient of the Hamiltonian with respect to the positions of the N monomers in the polymer
calculate_Hamiltonian_and_Gradient <- function(positions, types) {
  N <- nrow(positions)
  d <- ncol(positions)
  
  # Calculate the Lennard-Jones potential
  LJ_potential <- 0
  LJ_gradient <- matrix(0, nrow = N, ncol = d)
  for (i in 1:(N - 2)) {
    for (j in (i + 2):N) {
      if (types[i] == 'A' & types[j] == 'A'){c = 1}
      if (types[i] == 'A' & types[j] == 'B'){c = -0.5}
      if (types[i] == 'B' & types[j] == 'A'){c = -0.5}
      if (types[i] == 'B' & types[j] == 'B'){c = 0.5}
      # Determine the scaling factor for the Lennard-Jones potential
      # c <- ifelse(types[i] == types[j], 1, ifelse(types[i] == "A" & types[j] == "B", -1/2, 1/2))
      
      # Calculate the distance between two monomers
      r <- sqrt(sum((positions[i, ] - positions[j, ])^2))
      
      # Update the Lennard-Jones potential and its gradient
      LJ_potential <- LJ_potential + calculate_LJ_potential(r, 1, c)
      LJ_gradient[j, ] <- LJ_gradient[j, ] - (positions[i, ] - positions[j, ]) * 12 * (1 / r)^14 * (1 - c / 2 * (1 / r)^6)
    }
  }
  
  # Calculate the bond angle potential
  bond_potential <- 0
  bond_gradient <- matrix(0, nrow = N, ncol = d)
  for (i in 1:(N - 2)) {
    # Calculate the bond angle
    bond_vector_1 <- positions[i + 1, ] - positions[i, ]
    bond_vector_2 <- positions[i + 2, ] - positions[i + 1, ]
    bond_angle <- acos(sum(bond_vector_1 * bond_vector_2) / (sqrt(sum(bond_vector_1^2)) * sqrt(sum(bond_vector_2^2))))
    
    # Update the bond angle potential and its gradient
    bond_potential <- bond_potential + calculate_bond_potential(bond_angle)
    bond_gradient[i + 1, ] <- bond_gradient[i + 1, ] + (bond_vector_2 / sqrt(sum(bond_vector_2^2)) - bond_vector_1 * sum(bond_vector_1 * bond_vector_2) / (sqrt(sum(bond_vector_1^2)) * sum(bond_vector_2^2))) * (1 / 4 * sin(bond_angle))
  }
  
  # Calculate the Hamiltonian and its gradient
  Hamiltonian <- LJ_potential + bond_potential
  gradient <- LJ_gradient + bond_gradient
  
  list(Hamiltonian = Hamiltonian, gradient = gradient)
}

                                                                                                      