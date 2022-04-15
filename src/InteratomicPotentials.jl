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
include("NeighborList/nnlist.jl")
include("Parameters/param.jl")
include("PotentialTypes/types.jl")

# Export energies, forces, virial evaluations 
export energy_and_force, potential_energy, force, virial, virial_stress

export get_parameters, set_parameters, serialize_parameters, deserialize_parameters
export get_hyperparameters, set_hyperparameters, serialize_hyperparameters, deserialize_hyperparameters

end
