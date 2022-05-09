################################################################################
# Types of Empirical Potentials
################################################################################

include("ACE/ace.jl")
include("SNAP/snap.jl")

export ACE, ACEParams, SNAP, SNAPParams

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::BasisPotential)
    B, dB, _ = evaluate_full(A, p.basis_params)
    e = sum(B) ⋅ p.coefficients * ENERGY_UNIT
    f = [SVector{3}(di' * p.coefficients) for di in dB] * FORCE_UNIT
    (; e, f)
end

function virial_stress(A::AbstractSystem, p::BasisPotential)
    sum(SVector{6}(di ⋅ p.coefficients) for di in evaluate_basis_v(A, p.basis_params)) * ENERGY_UNIT
end

