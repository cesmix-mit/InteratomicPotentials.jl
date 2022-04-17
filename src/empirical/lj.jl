############################## Lennard Jones ###################################
struct LennardJones{T<:AbstractFloat} <: EmpiricalPotential{Parameter{(:ϵ, :σ)},Parameter{(:rcutoff,)}}
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

potential_energy(R::AbstractFloat, lj::LennardJones) = 4lj.ϵ * ((lj.σ / R)^12 - (lj.σ / R)^6)
force(R::AbstractFloat, r::SVector{3}, lj::LennardJones) = (24lj.ϵ * (2(lj.σ / R)^12 - (lj.σ / R)^6) / R^2)r
