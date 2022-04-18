"""
    LennardJones{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:ϵ, :σ)},NamedTuple{(:rcutoff,)}}

# TODO (Dallas): description of LennardJones
"""
struct LennardJones{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:ϵ, :σ)},NamedTuple{(:rcutoff,)}}
    ϵ::T
    σ::T
    rcutoff::T
    species::Tuple
    function LennardJones(ϵ, σ, rcutoff, species)
        ϵ, σ, rcutoff = promote(austrip(ϵ), austrip(σ), austrip(rcutoff))
        new{typeof(rcutoff)}(ϵ, σ, rcutoff, Tuple(species))
    end
    function LennardJones{T}(; ϵ::T, σ::T, rcutoff::T, species::Tuple) where {T}
        new{T}(ϵ, σ, rcutoff, species)
    end
end

get_rcutoff(lj::LennardJones) = lj.rcutoff
get_species(lj::LennardJones) = lj.species

potential_energy(R::AbstractFloat, lj::LennardJones) = 4lj.ϵ * ((lj.σ / R)^12 - (lj.σ / R)^6)
force(R::AbstractFloat, r::SVector{3}, lj::LennardJones) = (24lj.ϵ * (2(lj.σ / R)^12 - (lj.σ / R)^6) / R^2)r
