module InteratomicBasisPotentials

using InteratomicPotentials
using StaticArrays 
using LinearAlgebra
using AtomsBase
using Unitful
using UnitfulAtomic

# Include types
include("basis_potentials.jl")

#include SNAP
include("SNAP/snap.jl")
export SNAP, SNAPParams, get_num_snap_coeffs  # Export SNAP

include("ACE/ace.jl")
export RPI, RPIParams, get_rpi # ACE

# Export Basis set evaluations 
export evaluate_basis, evaluate_basis_d, evaluate_basis_v, evaluate_full
export energy_and_force, virial_stress
end # module
