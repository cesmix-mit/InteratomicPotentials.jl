##############################   Coulomb  ###################################
mutable struct Coulomb <: EmpiricalPotential
    q_1::Float64
    q_2::Float64
    ϵ0::Float64
end

function Coulomb()
    #ToDO
    return Coulomb(1.0, 1.0, 1.0)
end

function get_trainable_params(c::Coulomb)
    p = Parameter{}(())
    return p
end

function get_nontrainable_params(c::Coulomb)
    return (q_1 = c.q_1, q_2 = c.q_2, ϵ0 = c.ϵ0)
end

############################# Energies ##########################################

function potential_energy(r::Position, p::Coulomb)
    return p.q_1 * p.q_2 / (4.0 * π * p.ϵ0 * norm(r))
end

############################### Forces ##########################################

function force(r:: Position, p::Coulomb)
    d = norm(r)
    return p.q_1 * p.q_2 / (4.0 * π * p.ϵ0 * d^2) .* [r.x, r.y, r.z] ./ d
end

##############################   Gradients  ###################################

function grad_potential_energy(r::Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dpdnull = 0.0)
end
function grad_force(r:: Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dfdnull =  0 .* [r.x, r.y, r.z] ./ d )
end

function grad_virial(r::Position, p::Coulomb)
    println("The Coulomb potential has no trainable parameters")
    return (dvdnull =  0 .* [r.x, r.y, r.z] ./ d )
end

