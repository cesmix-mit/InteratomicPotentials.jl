##############################   Born-Mayer  ###################################
struct BornMayer <: EmpiricalPotential
    A
    ρ
    rcutoff
    species::AbstractVector
end
function BornMayer(A::Unitful.Energy, ρ::Unitful.Length, rcutoff::Unitful.Length, species::AbstractVector)
    BornMayer(austrip(A), austrip(ρ), austrip(rcutoff), species)
end

get_parameters(bm::BornMayer) = Parameter{(:A, :ρ)}((bm.A, bm.ρ))
set_parameters(p::Parameter{(:A, :ρ)}, bm::BornMayer) = BornMayer(p.A, p.ρ, bm.rcutoff, bm.species)

deserialize_parameters(p::Parameter{(:A, :ρ)}, bm::BornMayer) = [p.A, p.ρ]
serialize_parameters(p::Vector, bm::BornMayer) = Parameter{(:A, :ρ)}((p[1], p[2]))

get_hyperparameters(bm::BornMayer) = Parameter{(:rcutoff,)}((bm.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, bm::BornMayer) = BornMayer(bm.A, bm.ρ, p.rcutoff, bm.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, bm::BornMayer) = [p.rcutoff]
serialize_hyperparameters(p::Vector, bm::BornMayer) = Parameter{(:rcutoff,)}((p[1],))

# ##############################   Energy  ###################################

function potential_energy(R::AbstractFloat, bm::BornMayer)
    bm.A * exp(-R / bm.ρ)
end


# ##############################   Force   ###################################

function force(R::AbstractFloat, r::SVector{3}, bm::BornMayer)
    (bm.A * exp(-R / bm.ρ) / (bm.ρ * R))r
end

# ##############################   Gradients  ###################################

function grad_potential_energy(R::AbstractFloat, bm::BornMayer)
    (dpdA=exp(-R / bm.ρ),
        dpdρ=bm.A * R * exp(-R / bm.ρ) / bm.ρ^2)
end

function grad_force(R::AbstractFloat, r::SVector{3}, bm::BornMayer)
    (dfdA=(bm.ρ * exp(-R / bm.ρ) / (bm.ρ * R))r,
        dfdρ=(bm.A * (R - bm.ρ) * exp(-R / bm.ρ) / (bm.ρ^3 * R))r)
end

function grad_virial(r::SVector{3}, bm::BornMayer)
    dfdA, dfdρ = grad_force(r, bm)
    (dvdA=dfdA ⋅ r, dvdρ=dfdρ ⋅ r)
end
