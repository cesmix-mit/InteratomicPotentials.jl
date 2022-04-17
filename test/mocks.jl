struct MockArbitraryPotential <: ArbitraryPotential
    seed::Float64
    species::Tuple
end
function InteratomicPotentials.energy_and_force(A::AbstractSystem, p::MockArbitraryPotential)
    (
        e=(p.seed)u"hartree",
        f=fill((@SVector[p.seed + 1, p.seed + 2, p.seed + 3])u"hartree/bohr", length(A))
    )
end
function InteratomicPotentials.virial_stress(A::AbstractSystem, p::MockArbitraryPotential)
    (@SVector [p.seed - 1, p.seed - 2, p.seed - 3, p.seed - 4, p.seed - 5, p.seed - 6])u"hartree"
end

struct MockEmpiricalPotential <: EmpiricalPotential{Parameter{},Parameter{}}
    seed::Float64
    rcutoff::Float64
    species::Tuple
end
function InteratomicPotentials.potential_energy(R::AbstractFloat, p::MockEmpiricalPotential)
    R^2 * p.seed
end
function InteratomicPotentials.force(R::AbstractFloat, r::SVector{3}, p::MockEmpiricalPotential)
    R^2 * p.seed * r
end
