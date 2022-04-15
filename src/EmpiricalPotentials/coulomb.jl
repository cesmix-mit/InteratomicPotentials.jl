##############################   Coulomb  ###################################
struct Coulomb{T<:AbstractFloat} <: EmpiricalPotential
    q₁::T
    q₂::T
    rcutoff::T
    species::Vector{Symbol}
end
function Coulomb(q₁::Unitful.Charge, q₂::Unitful.Charge, rcutoff::Unitful.Length, species::AbstractVector{Symbol})
    Coulomb(austrip(q₁), austrip(q₂), austrip(rcutoff), collect(species))
end

get_parameters(::Coulomb) = Parameter{}()
set_parameters(c::Coulomb, ::Parameter{}) = copy(c)

serialize_parameters(c::Coulomb) = collect(get_parameters(c))
deserialize_parameters(c::Coulomb, p::AbstractVector) = set_parameters(c, Parameter{}(p))

get_hyperparameters(c::Coulomb) = Parameter{:rcutoff}((c.rcutoff,))
set_hyperparameters(c::Coulomb, p::Parameter{(:rcutoff,)}) = Coulomb(c.q₁, c.q₂, p.rcutoff, c.species)

serialize_hyperparameters(c::Coulomb) = collect(get_hyperparameters(c))
deserialize_hyperparameters(c::Coulomb, p::AbstractVector) = set_hyperparameters(c, Parameter{(:rcutoff,)}(p))

potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q₁ * c.q₂ / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q₁ * c.q₂ / R^3)r
