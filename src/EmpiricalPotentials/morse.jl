############################## Morse ###################################
struct Morse{T<:AbstractFloat} <: EmpiricalPotential
    D::T
    α::T
    σ::T
    rcutoff::T
    species::Vector{Symbol}
end
function Morse(D::Unitful.Energy, α::AbstractFloat, σ::Unitful.Length, rcutoff::Unitful.Length, species::AbstractVector{Symbol})
    Morse(austrip(D), α, austrip(σ), austrip(rcutoff), collect(species))
end

get_parameters(m::Morse) = Parameter{(:D, :α, :σ)}((m.D, m.α, m.σ))
set_parameters(m::Morse, p::Parameter{(:D, :α, :σ)}) = Morse(p.D, p.α, p.σ, m.rcutoff, m.species)

serialize_parameters(m::Morse) = collect(get_parameters(m))
deserialize_parameters(m::Morse, p::AbstractVector) = set_parameters(m, Parameter{(:D, :α, :σ)}(p))

get_hyperparameters(m::Morse) = Parameter{(:rcutoff,)}((m.rcutoff,))
set_hyperparameters(m::Morse, p::Parameter{(:rcutoff,)}) = Morse(m.D, m.α, m.σ, p.rcutoff, m.species)

serialize_hyperparameters(m::Morse) = collect(get_hyperparameters(m))
deserialize_hyperparameters(m::Morse, p::AbstractVector) = set_hyperparameters(m, Parameter{(:rcutoff,)}(p))

_morse_exp(R::AbstractFloat, m::Morse) = exp(-m.α * (R - m.σ))

potential_energy(R::AbstractFloat, m::Morse) = m.D * (1 - _morse_exp(R, m))^2
force(R::AbstractFloat, r::SVector{3}, m::Morse) = (2 * m.D * m.α * (1 - _morse_exp(R, m)) * _morse_exp(R, m) / R)r
