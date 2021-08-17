############################## Lennard Jones ###################################
mutable struct LennardJones <: EmpiricalPotential
    ϵ::Real
    σ::Real
end

function LennardJones()
    #ToDO
    return LennardJones(1.0, 1.0)
end

function get_trainable_params(lj::LennardJones)
    p = Parameter{:ϵ, :σ}((lj.ϵ, lj.σ))
    return p
end

function get_nontrainable_params(lj::LennardJones)
    p = Parameter{}(())
    return p
end

############################# Energies ##########################################

function potential_energy(r::Position, p::LennardJones)
    d = p.σ / norm(r)
    return 4.0 * p.ϵ * ( d^12 - d^6 )
end

############################### Forces ##########################################

function force(r::Position, p::LennardJones)
    d = norm(r)
    return 24.0 * p.ϵ * ( 2.0 * ( p.σ / d )^12 -  ( p.σ / d)^6 ) .* [r.x, r.y, r.z] ./ d ./ d
end

