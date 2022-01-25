import Base.+
import Base./
import Base.*
struct LinearCombinationPotential <: MixedPotential
    potentials   :: Vector{ArbitraryPotential}
    coefficients :: Vector{Real}
end


+(p1::ArbitraryPotential, p2::ArbitraryPotential) = LinearCombinationPotential([p1, p2], [1.0, 1.0])
-(p1::ArbitraryPotential, p2::ArbitraryPotential) = LinearCombinationPotential([p1, p2], [1.0, -1.0])

+(mp::LinearCombinationPotential, p1::ArbitraryPotential) = LinearCombinationPotential(push!(mp.potentials, p1), push!(mp.coefficients, 1.0))
-(mp::LinearCombinationPotential, p1::ArbitraryPotential) = LinearCombinationPotential(push!(mp.potentials, p1), push!(mp.coefficients, -1.0))

+(mp1::LinearCombinationPotential, mp2::LinearCombinationPotential) = LinearCombinationPotential(append!(mp1.potentials, mp2.potentials), append!(mp1.coefficients, mp2.coefficients))
-(mp1::LinearCombinationPotential, mp2::LinearCombinationPotential) = LinearCombinationPotential(append!(mp1.potentials, mp2.potentials), append!(mp1.coefficients, -1.0*mp2.coefficients))

*(a :: Real, p::ArbitraryPotential) = LinearCombinationPotential([p], [a])
*(a :: Real, mp::LinearCombinationPotential) = LinearCombinationPotential(mp.potentials, a*mp.coefficients)
*(p1::ArbitraryPotential, p2::ArbitraryPotential) = LinearCombinationPotential([p1, p2], [1.0, -1.0])

/(p::ArbitraryPotential, a :: Real) = LinearCombinationPotential([p], [1.0/a])
/(mp::LinearCombinationPotential, a :: Real) = LinearCombinationPotential(mp.potentials, mp.coefficients/a)

function potential_energy(sys::AbstractSystem, mp::LinearCombinationPotential) 
    sum( [c*potential_energy(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients) ] )
end

function force(sys::AbstractSystem, mp::LinearCombinationPotential) 
    sum( [c.*force(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients) ] )
end

function virial(sys::AbstractSystem, mp::LinearCombinationPotential) 
    sum( [c*virial(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients) ] )
end

function virial_stress(sys::AbstractSystem, mp::LinearCombinationPotential) 
    sum( [c.*virial_stress(sys, p) for (p, c) in zip(mp.potentials, mp.coefficients) ] )
end
