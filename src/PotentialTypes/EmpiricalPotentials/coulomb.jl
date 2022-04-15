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

get_parameters(c::Coulomb) = Parameter{}(())
set_parameters(p::Parameter{}, c::Coulomb) = c

deserialize_parameters(p::Parameter{()}, c::Coulomb) = []
serialize_parameters(p::Vector, c::Coulomb) = Parameter{()}(())

get_hyperparameters(c::Coulomb) = Parameter{:rcutoff}((c.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = Coulomb(c.q₁, c.q₂, p.rcutoff, c.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = [p.rcutoff]
serialize_hyperparameters(p::Vector, c::Coulomb) = Parameter{(:rcutoff,)}((p[1],))

potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q₁ * c.q₂ / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q₁ * c.q₂ / R^3)r
