abstract type BasisPotential <: AbstractPotential end
abstract type BasisParameters end

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################
import InteratomicPotentials: energy_and_force, virial_stress, ENERGY_UNIT, FORCE_UNIT

function energy_and_force(A::AbstractSystem, p::BasisPotential)
    B, dB, _ = evaluate_full(A, p.basis_params)
    e = sum(B) ⋅ p.coefficients * ENERGY_UNIT
    f = [SVector{3}(di' * p.coefficients) for di in dB] * FORCE_UNIT
    (; e, f)
end

function virial_stress(A::AbstractSystem, p::BasisPotential)
    sum(SVector{6}(di ⋅ p.coefficients) for di in evaluate_basis_v(A, p.basis_params)) * ENERGY_UNIT
end
