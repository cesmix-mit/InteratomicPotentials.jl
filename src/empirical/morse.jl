"""
    Morse{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:D, :α, :σ)},NamedTuple{(:rcutoff,)}}

The Morse potential is a simple two-body intermolecular potential with three typical parameters, primarily describing the interaction of neutral atoms with more complex bond interactions. Formally, the interaction between two atoms at a distance, ``r``, is given by 

```math
\\begin{equation}
V_{M}(r; D, \\alpha, \\sigma, rcutoff, species) =
    \\begin{cases} 
    0 & r > rcutoff \\\\
    D \\left( 1 - e^{\\alpha(r - \\sigma)}\right)^2 & r < rcutoff.
    \\end{cases}
\\end{equation}
```

Users must supply three parameters, ``D`` (units of energy), ``\\alpha`` (units of inverse distance), ''\\sigma'' (units of distance), and radial cutoff (units of distance).
"""
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
