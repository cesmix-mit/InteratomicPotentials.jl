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

potential_energy(R::AbstractFloat, bm::BornMayer) = bm.A * exp(-R / bm.ρ)
force(R::AbstractFloat, r::SVector{3}, bm::BornMayer) = (bm.A * exp(-R / bm.ρ) / (bm.ρ * R))r
