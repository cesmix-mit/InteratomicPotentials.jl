############################## Lennard Jones ###################################
struct LennardJones{T<:AbstractFloat} <: EmpiricalPotential
    ϵ::T
    σ::T
    rcutoff::T
    species::Tuple
    function LennardJones(ϵ, σ, rcutoff, species)
        ϵ, σ, rcutoff = promote(austrip(ϵ), austrip(σ), austrip(rcutoff))
        new{typeof(rcutoff)}(ϵ, σ, rcutoff, Tuple(species))
    end
end

get_rcutoff(lj::LennardJones) = lj.rcutoff
get_species(lj::LennardJones) = lj.species

get_parameters(lj::LennardJones) = Parameter{(:ϵ, :σ)}((lj.ϵ, lj.σ))
set_parameters(lj::LennardJones, p::Parameter{(:ϵ, :σ)}) = LennardJones(p.ϵ, p.σ, lj.rcutoff, lj.species)

serialize_parameters(lj::LennardJones) = collect(get_parameters(lj))
deserialize_parameters(lj::LennardJones, p::AbstractVector) = set_parameters(lj, Parameter{(:ϵ, :σ)}(p))

get_hyperparameters(lj::LennardJones) = Parameter{(:rcutoff,)}((lj.rcutoff,))
set_hyperparameters(lj::LennardJones, p::Parameter{(:rcutoff,)}) = LennardJones(lj.ϵ, lj.σ, p.rcutoff, lj.species)

serialize_hyperparameters(lj::LennardJones) = collect(get_hyperparameters(lj))
deserialize_hyperparameters(lj::LennardJones, p::AbstractVector) = set_hyperparameters(lj, Parameter{(:rcutoff,)}(p))

potential_energy(R::AbstractFloat, lj::LennardJones) = 4lj.ϵ * ((lj.σ / R)^12 - (lj.σ / R)^6)
force(R::AbstractFloat, r::SVector{3}, lj::LennardJones) = (24lj.ϵ * (2(lj.σ / R)^12 - (lj.σ / R)^6) / R^2)r
