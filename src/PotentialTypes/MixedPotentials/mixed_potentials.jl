struct LinearCombinationPotential <: MixedPotential
    potentials::Vector{ArbitraryPotential}
    coefficients::Vector{AbstractFloat}
end
LinearCombinationPotential(p::ArbitraryPotential) = LinearCombinationPotential([p], [1])
LinearCombinationPotential(p::LinearCombinationPotential) = p

Base.:+(p1::ArbitraryPotential, p2::ArbitraryPotential) = LinearCombinationPotential(p1) + LinearCombinationPotential(p2)
Base.:+(p::ArbitraryPotential) = +1 * p
Base.:-(p1::ArbitraryPotential, p2::ArbitraryPotential) = LinearCombinationPotential(p1) + -LinearCombinationPotential(p2)
Base.:-(p::ArbitraryPotential) = -1 * p
Base.:*(a::Real, p::ArbitraryPotential) = a * LinearCombinationPotential(p)
Base.:/(p::ArbitraryPotential, a::Real) = LinearCombinationPotential(p) / a

Base.:+(mp1::LinearCombinationPotential, mp2::LinearCombinationPotential) = LinearCombinationPotential(append!(mp1.potentials, mp2.potentials), append!(mp1.coefficients, mp2.coefficients))
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
