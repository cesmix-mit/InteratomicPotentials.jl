############################## Lennard Jones ###################################
struct LennardJones <: EmpiricalPotential
    ϵ::AbstractFloat
    σ::AbstractFloat
    rcutoff::AbstractFloat
end

get_trainable_params(lj::LennardJones) = Parameter{:ϵ,:σ}((lj.ϵ, lj.σ))

get_nontrainable_params(lj::LennardJones) = Parameter{}(())

############################# Energies ##########################################

function potential_energy(R::AbstractFloat, p::LennardJones)
    d = p.σ / R
    4.0 * p.ϵ * (d^12 - d^6)
end

############################### Forces ##########################################

force(r::SVector{3,<:AbstractFloat}, p::LennardJones) = force(r, norm(r), p)

function force(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, p::LennardJones)
    SVector(24.0 * p.ϵ * (2.0 * (p.σ / R)^12 - (p.σ / R)^6) .* r ./ R ./ R)
end
############################## Gradients ########################################
function grad_potential_energy(r::SVector{3,<:AbstractFloat}, p::LennardJones)
    d = p.σ / norm(r)
    (dpdϵ = 4.0 * (d^12 - d^6), dpdσ = 24.0 * p.ϵ / p.σ * (2 * d^12 - d^6))
end

function grad_force(r::SVector{3,<:Real}, p::LennardJones)
    d = norm(r)
    (dfdϵ = 48.0 * ((p.σ / d)^13 - 0.5 * (p.σ / d)^7) .* r ./ d,
        dfdσ = 144.0 * p.ϵ / p.σ * (4.0 * p.σ^12 / d^13 - p.σ^6 / d^7) .* r ./ d)
end

function grad_virial(r::SVector{<:Real}, p::LennardJones)
    df = grad_force(r, p)
    dfde = df[1]
    dfdsig = df[2]

    (dvdϵ = dfde[1] * r[1] + dfde[2] * r[2] + dfde[3] * r[3],
        dvdσ = dfdsig[1] * r[1] + dfdsig[2] * r[2] + dfdsig[3] * r[3])
end
