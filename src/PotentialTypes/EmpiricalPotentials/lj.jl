############################## Lennard Jones ###################################
struct LennardJones <: EmpiricalPotential
    ϵ
    σ
    rcutoff
    species::AbstractVector
end

get_parameters(lj::LennardJones) = Parameter{(:ϵ, :σ)}((lj.ϵ, lj.σ))
set_parameters(p::Parameter{(:ϵ, :σ)}, lj::LennardJones) = LennardJones(p.ϵ, p.σ, lj.rcutoff, lj.species)

deserialize_parameters(p::Parameter{(:ϵ, :σ)}, lj::LennardJones) = [p.ϵ, p.σ]
serialize_parameters(p::Vector, lj::LennardJones) = Parameter{(:ϵ, :σ)}((p[1], p[2]))

get_hyperparameters(lj::LennardJones) = Parameter{(:rcutoff,)}((lj.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, lj::LennardJones) = LennardJones(lj.ϵ, lj.σ, p.rcutoff, lj.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, lj::LennardJones) = [p.rcutoff]
serialize_hyperparameters(p::Vector, lj::LennardJones) = Parameter{(:rcutoff,)}((p[1],))

############################# Energies ##########################################

function potential_energy(R::AbstractFloat, p::LennardJones)
    d = p.σ / R
    4.0 * p.ϵ * (d^12 - d^6)
end

############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3}, p::LennardJones)
    SVector(24.0 * p.ϵ * (2.0 * (p.σ / R)^12 - (p.σ / R)^6) .* r ./ R ./ R)
end

############################## Gradients ########################################
function grad_potential_energy(r::SVector{3}, p::LennardJones)
    d = p.σ / norm(r)
    (dpdϵ = 4.0 * (d^12 - d^6),
        dpdσ = 24.0 * p.ϵ / p.σ * (2 * d^12 - d^6))
end

function grad_force(r::SVector{3}, p::LennardJones)
    d = norm(r)
    (dfdϵ = 48.0 * ((p.σ / d)^13 - 0.5 * (p.σ / d)^7) .* r ./ d,
        dfdσ = 144.0 * p.ϵ / p.σ * (4.0 * p.σ^12 / d^13 - p.σ^6 / d^7) .* r ./ d)
end

function grad_virial(r::SVector{3}, p::LennardJones)
    dfdϵ, dfdσ = grad_force(r, p)
    (dvdϵ = dfdϵ ⋅ r, dvdσ = dfdσ ⋅ r)
end
