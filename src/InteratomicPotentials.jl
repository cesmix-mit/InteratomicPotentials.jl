module InteratomicPotentials

using AtomsBase
using LinearAlgebra
using StaticArrays
using Unitful
using UnitfulAtomic

using Distances
using NearestNeighbors

# Import files
include("unit_convention.jl")
include("constants.jl")
include("nnlist.jl")
include("types.jl")

# Export energies, forces, virial evaluations 
export energy_and_force, potential_energy, force, virial, virial_stress

# Export parameter manipulation functions
export get_parameters, set_parameters, serialize_parameters, deserialize_parameters
export get_hyperparameters, set_hyperparameters, serialize_hyperparameters, deserialize_hyperparameters

end
