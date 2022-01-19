##############################   Coulomb  ###################################
mutable struct Coulomb <: EmpiricalPotential
    q_1::AbstractFloat
    q_2::AbstractFloat
    ϵ0::AbstractFloat
    rcutoff::AbstractFloat
end

get_trainable_params(c::Coulomb) = Parameter{}(())

get_nontrainable_params(c::Coulomb) = Parameter{:q1, :q2, :ϵ0, :rcutoff}((c.q_1, c.q_2, c.ϵ0, c.rcutoff))

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
