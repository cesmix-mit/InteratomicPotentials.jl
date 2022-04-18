############################## Morse ###################################
struct Morse{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:D, :α, :σ)},NamedTuple{(:rcutoff,)}}
    D::T
    α::T
    σ::T
    rcutoff::T
    species::Tuple
    function Morse(D, α, σ, rcutoff, species)
        D, α, σ, rcutoff = promote(austrip(D), α, austrip(σ), austrip(rcutoff))
        new{typeof(rcutoff)}(D, α, σ, rcutoff, Tuple(species))
    end
    function Morse{T}(; D::T, α::T, σ::T, rcutoff::T, species::Tuple) where {T}
        new{T}(D, α, σ, rcutoff, species)
    end
end

get_rcutoff(m::Morse) = m.rcutoff
get_species(m::Morse) = m.species

_morse_exp(R::AbstractFloat, m::Morse) = exp(-m.α * (R - m.σ))

potential_energy(R::AbstractFloat, m::Morse) = m.D * (1 - _morse_exp(R, m))^2
force(R::AbstractFloat, r::SVector{3}, m::Morse) = (2 * m.D * m.α * (1 - _morse_exp(R, m)) * _morse_exp(R, m) / R)r
