struct LinearCombinationPotential <: MixedPotential
    potentials::Vector{ElementalPotential}
    coefficients::Vector{AbstractFloat}
end
LinearCombinationPotential(p::ElementalPotential) = LinearCombinationPotential([p], [1])
LinearCombinationPotential(p::LinearCombinationPotential) = p

Base.:+(p::ElementalPotential) = +1 * p
Base.:+(mp::LinearCombinationPotential) = mp
Base.:+(p1::ElementalPotential, p2::ElementalPotential) = LinearCombinationPotential(p1) + LinearCombinationPotential(p2)

Base.:-(p::ElementalPotential) = -1 * p
Base.:-(mp::LinearCombinationPotential) = LinearCombinationPotential(mp.potentials, -1.0*mp.coefficients)
Base.:-(p1::ElementalPotential, p2::ElementalPotential) = LinearCombinationPotential(p1) + -LinearCombinationPotential(p2)

Base.:*(a::Real, p::ElementalPotential) = a * LinearCombinationPotential(p)
Base.:/(p::ElementalPotential, a::Real) = LinearCombinationPotential(p) / a

Base.:+(mp1::LinearCombinationPotential, mp2::LinearCombinationPotential) = LinearCombinationPotential(append!(mp1.potentials, mp2.potentials), append!(mp1.coefficients, mp2.coefficients))
Base.:-(mp1::LinearCombinationPotential, mp2::LinearCombinationPotential) = mp1 + -mp2
Base.:*(a::Real, mp::LinearCombinationPotential) = LinearCombinationPotential(mp.potentials, a * mp.coefficients)
Base.:/(mp::LinearCombinationPotential, a::Real) = LinearCombinationPotential(mp.potentials, mp.coefficients / a)

function potential_energy(sys::AbstractSystem, mp::LinearCombinationPotential)
    sum(c * potential_energy(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients))
end

function force(sys::AbstractSystem, mp::LinearCombinationPotential)
    sum(c .* force(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients))
end

function virial(sys::AbstractSystem, mp::LinearCombinationPotential)
    sum(c * virial(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients))
end

function virial_stress(sys::AbstractSystem, mp::LinearCombinationPotential)
    sum(c .* virial_stress(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients))
end
