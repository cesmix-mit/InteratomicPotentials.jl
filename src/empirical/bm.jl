##############################   Born-Mayer  ###################################
struct BornMayer{T<:AbstractFloat} <: EmpiricalPotential
    A::T
    ρ::T
    rcutoff::T
    species::Tuple
    function BornMayer(A, ρ, rcutoff, species)
        A, ρ, rcutoff = promote(austrip(A), austrip(ρ), austrip(rcutoff))
        new{typeof(rcutoff)}(A, ρ, rcutoff, Tuple(species))
    end
end

get_rcutoff(bm::BornMayer) = bm.rcutoff
get_species(bm::BornMayer) = bm.species

get_parameters(bm::BornMayer) = Parameter{(:A, :ρ)}((bm.A, bm.ρ))
set_parameters(bm::BornMayer, p::Parameter{(:A, :ρ)}) = BornMayer(p.A, p.ρ, bm.rcutoff, bm.species)

serialize_parameters(bm::BornMayer) = collect(get_parameters(bm))
deserialize_parameters(bm::BornMayer, p::AbstractVector) = set_parameters(bm, Parameter{(:A, :ρ)}(p))

get_hyperparameters(bm::BornMayer) = Parameter{(:rcutoff,)}((bm.rcutoff,))
set_hyperparameters(bm::BornMayer, p::Parameter{(:rcutoff,)}) = BornMayer(bm.A, bm.ρ, p.rcutoff, bm.species)

serialize_hyperparameters(bm::BornMayer) = collect(get_hyperparameters(bm))
deserialize_hyperparameters(bm::BornMayer, p::AbstractVector) = set_hyperparameters(bm, Parameter{(:rcutoff,)}(p))

potential_energy(R::AbstractFloat, bm::BornMayer) = bm.A * exp(-R / bm.ρ)
force(R::AbstractFloat, r::SVector{3}, bm::BornMayer) = (bm.A * exp(-R / bm.ρ) / (bm.ρ * R))r
