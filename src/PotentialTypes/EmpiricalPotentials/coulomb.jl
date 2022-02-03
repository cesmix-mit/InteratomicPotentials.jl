##############################   Coulomb  ###################################
mutable struct Coulomb <: EmpiricalPotential
    q_1
    q_2
    ϵ0
    rcutoff
    species :: AbstractVector
end

get_parameters(c::Coulomb) = Parameter{}(())
set_parameters(p::Parameter{}, c::Coulomb) = c

deserialize_parameters(p::Parameter{()}, c::Coulomb) = []
serialize_parameters(p::Vector, c::Coulomb) = Parameter{()}( () )

get_hyperparameters(c::Coulomb) = Parameter{:rcutoff}((c.rcutoff, ))
set_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = Coulomb(c.q_1, c.q_2, c.ϵ0, p.rcutoff, c.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, c::Coulomb) = [p.rcutoff]
serialize_hyperparameters(p::Vector, c::Coulomb) = Parameter{(:rcutoff, )}( (p[1],) )

# ############################# Energies ##########################################

function potential_energy(R::AbstractFloat, c::Coulomb)
    return c.q_1 * c.q_2 / (4.0 * π * c.ϵ0 * R)
end

# ############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, c::Coulomb)
    SVector(c.q_1 * c.q_2 / (4.0 * π * c.ϵ0 * R^2) .* r ./ R)
end

# ##############################   Gradients  ###################################

function grad_potential_energy(r::Vector{<:Real}, p::Coulomb)
    error("The Coulomb potential has no trainable parameters")
end
function grad_force(r:: Vector{<:Real}, p::Coulomb)
    error("The Coulomb potential has no trainable parameters")
end

function grad_virial(r::Vector{<:Real}, p::Coulomb)
    error("The Coulomb potential has no trainable parameters")
end
