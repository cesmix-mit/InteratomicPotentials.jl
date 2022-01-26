# Include types
abstract type BasisParameters end

#include SNAP
include("SNAP/snap.jl")

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################

function force(A::AbstractSystem, p::BasisPotential)
    d = evaluate_basis_d(A, p.basis_params)
    [SVector{3}(di' * p.coefficients) for di in d]
end

function potential_energy(A::AbstractSystem, p::BasisPotential)
    dot(evaluate_basis(A, p.basis_params), p.coefficients)
end

function virial_stress(A::AbstractSystem, p::BasisPotential)
    d = evaluate_basis_v(A, p.basis_params)
    SVector{6}(sum(dot(di, p.coefficients) for di in d))
end
