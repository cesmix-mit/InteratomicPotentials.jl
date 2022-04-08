############################## Lennard Jones ###################################
struct Morse <: EmpiricalPotential
    D
    α
    σ
    rcutoff
    species::AbstractVector
end
function Morse(D::Unitful.Energy, α::Real, σ::Unitful.Length, rcutoff::Unitful.Length, species::AbstractVector)
    Morse(austrip(D), α, austrip(σ), austrip(rcutoff), species)
end

get_parameters(m::Morse) = Parameter{(:D, :α, :σ)}((m.D, m.α, m.σ))
set_parameters(p::Parameter{(:D, :α, :σ)}, m::Morse) = Morse(p.D, p.α, p.σ, m.rcutoff, m.species)

deserialize_parameters(p::Parameter{(:D, :α, :σ)}, m::Morse) = [p.ϵ, p.σ]
serialize_parameters(p::Vector, m::Morse) = Parameter{(:ϵ, :σ)}((p[1], p[2]))

get_hyperparameters(m::Morse) = Parameter{(:rcutoff,)}((m.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, m::Morse) = Morse(m.D, m.α, m.σ, p.rcutoff, m.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, m::Morse) = [p.rcutoff]
serialize_hyperparameters(p::Vector, m::Morse) = Parameter{(:rcutoff,)}((p[1],))

############################# Energies ##########################################

function potential_energy(R::AbstractFloat, p::Morse)
    p.D * (1 - exp(-p.α * (R - p.σ)))^2
end

############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3}, p::Morse)
    e = exp(-p.α * (R - p.σ))
    (2 * p.D * p.α * (1 - e) * e / R)r
end

###############################   Gradients  ###################################

# TODO
