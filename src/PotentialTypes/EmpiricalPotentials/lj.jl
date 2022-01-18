############################## Lennard Jones ###################################
mutable struct LennardJones <: EmpiricalPotential
    ϵ::AbstractFloat
    σ::AbstractFloat
    rcutoff::AbstractFloat
end

function get_trainable_params(lj::LennardJones)
    p = Parameter{:ϵ,:σ}((lj.ϵ, lj.σ))
    return p
end

function get_nontrainable_params(lj::LennardJones)
    p = Parameter{}(())
    return p
end

############################# Energies ##########################################

function potential_energy(R::AbstractFloat, p::LennardJones)
    d = p.σ / R
    return 4.0 * p.ϵ * (d^12 - d^6)
end

############################### Forces ##########################################

function force(r::SVector{3,<:AbstractFloat}, p::LennardJones)
    d = norm(r)
    return SVector(24.0 * p.ϵ * (2.0 * (p.σ / d)^12 - (p.σ / d)^6) .* r ./ d ./ d)
end

function force(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, p::LennardJones)
    return SVector(24.0 * p.ϵ * (2.0 * (p.σ / R)^12 - (p.σ / R)^6) .* r ./ R ./ R)
end
############################## Gradients ########################################
function grad_potential_energy(r::SVector{3,<:AbstractFloat}, p::LennardJones)
    d = p.σ / norm(r)
    g = (dpdϵ = 4.0 * (d^12 - d^6), dpdσ = 24.0 * p.ϵ / p.σ * (2 * d^12 - d^6))

    return g
end

function grad_force(r::SVector{3,<:Real}, p::LennardJones)
    d = norm(r)
    g = (dfdϵ = 48.0 * ((p.σ / d)^13 - 0.5 * (p.σ / d)^7) .* r ./ d,
        dfdσ = 144.0 * p.ϵ / p.σ * (4.0 * p.σ^12 / d^13 - p.σ^6 / d^7) .* r ./ d)
    return g
end

function grad_virial(r::SVector{<:Real}, p::LennardJones)
    df = grad_force(r, p)
    dfde = df[1]
    dfdsig = df[2]

    v = (dvdϵ = dfde[1] * r[1] + dfde[2] * r[2] + dfde[3] * r[3],
        dvdσ = dfdsig[1] * r[1] + dfdsig[2] * r[2] + dfdsig[3] * r[3])
    return v

end
