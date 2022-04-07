##############################   Coulomb  ###################################
struct Coulomb <: EmpiricalPotential
    q_1
    q_2
    ϵ0
    rcutoff
    species::AbstractVector
end

get_parameters(c::Coulomb) = Parameter{}(())
set_parameters(p::Parameter{}, c::Coulomb) = c

deserialize_parameters(p::Parameter{()}, c::Coulomb) = []
serialize_parameters(p::Vector, c::Coulomb) = Parameter{()}(())

get_hyperparameters(c::Coulomb) = Parameter{:rcutoff}((c.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = Coulomb(c.q_1, c.q_2, c.ϵ0, p.rcutoff, c.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = [p.rcutoff]
serialize_hyperparameters(p::Vector, c::Coulomb) = Parameter{(:rcutoff,)}((p[1],))

# ############################# Energies ##########################################

function potential_energy(R::AbstractFloat, c::Coulomb)
    c.q_1 * c.q_2 / (4π * c.ϵ0 * R)
end

# ############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3}, c::Coulomb)
    (c.q_1 * c.q_2 / (4π * c.ϵ0 * R^3))r
end

# ##############################   Gradients  ###################################

grad_potential_energy(R::AbstractFloat, p::Coulomb) = error("The Coulomb potential has no trainable parameters")
grad_force(R::AbstractFloat, r::SVector{3}, p::Coulomb) = error("The Coulomb potential has no trainable parameters")
grad_virial(r::SVector{3}, p::Coulomb) = error("The Coulomb potential has no trainable parameters")
