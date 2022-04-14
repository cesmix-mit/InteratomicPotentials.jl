##############################   Coulomb  ###################################
struct Coulomb <: EmpiricalPotential
    q_1
    q_2
    rcutoff
    species::AbstractVector
end
function Coulomb(q_1::Unitful.Charge, q_2::Unitful.Charge, rcutoff::Unitful.Length, species::AbstractVector)
    Coulomb(austrip(q_1), austrip(q_2), austrip(rcutoff), species)
end

get_parameters(c::Coulomb) = Parameter{}(())
set_parameters(p::Parameter{}, c::Coulomb) = c

deserialize_parameters(p::Parameter{()}, c::Coulomb) = []
serialize_parameters(p::Vector, c::Coulomb) = Parameter{()}(())

get_hyperparameters(c::Coulomb) = Parameter{:rcutoff}((c.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = Coulomb(c.q_1, c.q_2, p.rcutoff, c.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = [p.rcutoff]
serialize_hyperparameters(p::Vector, c::Coulomb) = Parameter{(:rcutoff,)}((p[1],))

potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q_1 * c.q_2 / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q_1 * c.q_2 / R^3)r
