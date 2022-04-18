##############################   Coulomb  ###################################
struct Coulomb{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{()},NamedTuple{(:rcutoff,)}}
    q₁::T
    q₂::T
    rcutoff::T
    species::Tuple
    function Coulomb(q₁, q₂, rcutoff, species)
        q₁, q₂, rcutoff = promote(austrip(q₁), austrip(q₂), austrip(rcutoff))
        new{typeof(rcutoff)}(q₁, q₂, rcutoff, Tuple(species))
    end
    function Coulomb{T}(; q₁::T, q₂::T, rcutoff::T, species::Tuple) where {T}
        new{T}(q₁, q₂, rcutoff, species)
    end
end

get_rcutoff(c::Coulomb) = c.rcutoff
get_species(c::Coulomb) = c.species
potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q₁ * c.q₂ / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q₁ * c.q₂ / R^3)r
