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

############################## Gradients ########################################
function grad_potential_energy(r::Position, p::LennardJones)
    d = p.σ / norm(r)
    g = (dpdϵ = 4.0 * (d^12 - d^6), dpdσ = 24.0 * p.ϵ / p.σ * (2*d^12 - d^6) )

    return g
end

function grad_force(r::Position, p::LennardJones)
    d = norm(r)
    g = (dfdϵ = 48.0 * ( (p.σ / d)^13 - 0.5*(p.σ / d)^7 ) .* [r.x, r.y, r.z] ./ d, 
    dfdσ = 144.0 * p.ϵ / p.σ * ( 4.0*p.σ^12 / d^13 - p.σ^6 / d^7 ) .* [r.x, r.y, r.z] ./ d )
    return g 
end

function grad_virial(r::Position, p::LennardJones)
    df = grad_force(r, p)
    dfde = df[1]
    dfdsig = df[2]
    
    v = (dvdϵ = dfde[1]*r.x + dfde[2]*r.y + dfde[3]*r.z, 
    dvdσ = dfdsig[1]*r.x + dfdsig[2]*r.y + dfdsig[3]*r.z)
    return v

end