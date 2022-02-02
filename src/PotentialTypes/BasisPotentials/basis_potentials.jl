# Include types
abstract type BasisParameters end

#include SNAP
include("SNAP/snap.jl")
export SNAP, SNAPParams, get_num_snap_coeffs  # Export SNAP

# include("ACE/ace.jl")
# export RPI, RPIParams # ACE

# Export Basis set evaluations 
export evaluate_basis, evaluate_basis_d, evaluate_basis_v, evaluate_full

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::BasisPotential)
    B, dB, _ = evaluate_full(A, p.basis_params)
    e = sum(B) ⋅ p.coefficients
    f = [SVector{3}(di' * p.coefficients) for di in dB]
    (; e, f)
end

function virial_stress(A::AbstractSystem, p::BasisPotential)
    sum(SVector{6}(di ⋅ p.coefficients) for di in evaluate_basis_v(A, p.basis_params))
end
