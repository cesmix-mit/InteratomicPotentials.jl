abstract type BasisPotential <: ArbitraryPotential end
abstract type BasisParameters end

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################
import InteratomicPotentials.energy_and_force
import InteratomicPotentials.virial_stress
function InteratomicPotentials.energy_and_force(A::AbstractSystem, p::BasisPotential)
    B, dB, _ = evaluate_full(A, p.basis_params)
    e = sum(B) ⋅ p.coefficients
    f = [SVector{3}(di' * p.coefficients) for di in dB]
    (; e, f)
end

function virial_stress(A::AbstractSystem, p::BasisPotential)
    sum(SVector{6}(di ⋅ p.coefficients) for di in evaluate_basis_v(A, p.basis_params))
end
