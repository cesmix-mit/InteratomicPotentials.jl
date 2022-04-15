##############################   Coulomb  ###################################
struct Coulomb{T<:AbstractFloat} <: EmpiricalPotential
    q₁::T
    q₂::T
    rcutoff::T
    species::Tuple
    function Coulomb(q₁, q₂, rcutoff, species)
        q₁, q₂, rcutoff = promote(austrip(q₁), austrip(q₂), austrip(rcutoff))
        new{typeof(rcutoff)}(q₁, q₂, rcutoff, Tuple(species))
    end
end

get_parameters(::Coulomb) = Parameter{}()
set_parameters(c::Coulomb, ::Parameter{}) = c

serialize_parameters(c::Coulomb) = collect(get_parameters(c))
deserialize_parameters(c::Coulomb, p::AbstractVector) = set_parameters(c, Parameter{}(p))

get_hyperparameters(c::Coulomb) = Parameter{(:rcutoff,)}((c.rcutoff,))
set_hyperparameters(c::Coulomb, p::Parameter{(:rcutoff,)}) = Coulomb(c.q₁, c.q₂, p.rcutoff, c.species)

serialize_hyperparameters(c::Coulomb) = collect(get_hyperparameters(c))
deserialize_hyperparameters(c::Coulomb, p::AbstractVector) = set_hyperparameters(c, Parameter{(:rcutoff,)}(p))

potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q₁ * c.q₂ / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q₁ * c.q₂ / R^3)r
