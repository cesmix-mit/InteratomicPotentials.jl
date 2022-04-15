##############################   Born-Mayer  ###################################
struct BornMayer{T<:AbstractFloat} <: EmpiricalPotential
    A::T
    ρ::T
    rcutoff::T
    species::Vector{Symbol}
end
function BornMayer(A::Unitful.Energy, ρ::Unitful.Length, rcutoff::Unitful.Length, species::AbstractVector{Symbol})
    BornMayer(austrip(A), austrip(ρ), austrip(rcutoff), collect(species))
end

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
