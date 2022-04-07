############################## Lennard Jones ###################################
struct LennardJones <: EmpiricalPotential
    ϵ
    σ
    rcutoff
    species::AbstractVector
end
function LennardJones(ϵ::Unitful.Energy, σ::Unitful.Length, rcutoff::Unitful.Length, species::AbstractVector)
    LennardJones(austrip(ϵ), austrip(σ), austrip(rcutoff), species)
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
    d⁶ = (p.σ / R)^6
    4p.ϵ * (d⁶^2 - d⁶)
end

############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3}, p::LennardJones)
    d⁶ = (p.σ / R)^6
    (24p.ϵ * (2d⁶^2 - d⁶) / R^2)r
end

############################## Gradients ########################################

function grad_potential_energy(r::SVector{3}, p::LennardJones)
    d⁶ = (p.σ / norm(r))^6
    (dpdϵ=4(d⁶^2 - d⁶),
        dpdσ=24p.ϵ * (2d⁶^2 - d⁶) / p.σ)
end

function grad_force(r::SVector{3}, p::LennardJones)
    R = norm(r)
    d⁶ = (p.σ / R)^6
    (dfdϵ=(24(2d⁶^2 - d⁶) / R^2)r,
        dfdσ=(144p.ϵ * (4d⁶^2 - d⁶) / (p.σ * R^2))r)
end

function grad_virial(r::SVector{3}, p::LennardJones)
    dfdϵ, dfdσ = grad_force(r, p)
    (dvdϵ=dfdϵ ⋅ r, dvdσ=dfdσ ⋅ r)
end
