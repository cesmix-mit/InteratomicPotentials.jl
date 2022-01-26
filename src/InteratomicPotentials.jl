module InteratomicPotentials

using StaticArrays
using LAMMPS
using LinearAlgebra
using AtomsBase
using Unitful
using UnitfulAtomic

# Import files
include("Utilities/utils.jl")
include("PotentialTypes/types.jl")

# Export Neighbor list tools
export NeighborList, neighborlist

# Export Potentials
export EmpiricalPotential, LennardJones, BornMayer, Coulomb, ZBL
export BasisPotential, BasisParameters
export SNAP, SNAPParams, get_num_snap_coeffs  # Export SNAP 

# Export energies, forces, virial evaluations 
export energy_and_force, potential_energy, force, virial, virial_stress
export grad_potential_energy, grad_force, grad_virial, grad_virial_stress

# Export Basis set evaluations 
export evaluate_basis, evaluate_basis_d, evaluate_basis_v, evaluate_full


end
