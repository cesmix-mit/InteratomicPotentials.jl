"""
    LennardJones{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:ϵ, :σ)},NamedTuple{(:rcutoff,)}}

The Lennard-Jones potential is a simple two-body intermolecular potential with two typical parameters, primarily describing the interaction of neutral atoms. Formally, the interaction between two atoms at a distance, ``r``, is given by 

```math
\begin{equation}
V_{LJ}(r; \\epsilon, \\sigma, rcutoff, species) =
    \begin{cases} 
        0 & r > rcutoff \\
        4\\epsilon \\lbrack \\frac{\\sigma}{r})^{12} - (\\frac{\\sigma}{r})^6 \\rbrack & r < rcutoff.
    \end{cases}
\end{equation}
```

Users must supply two parameters, ``\\epsilon`` (units of energy), ``\\sigma`` (units of distance), and radial cutoff (units of distance).
"""
struct LennardJones{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:ϵ, :σ)},NamedTuple{(:rcutoff,)}}
    ϵ::T
    σ::T
    rcutoff::T
    species::Tuple
    function LennardJones(ϵ, σ, rcutoff, species)
        ϵ, σ, rcutoff = promote(austrip(ϵ), austrip(σ), austrip(rcutoff))
        new{typeof(rcutoff)}(ϵ, σ, rcutoff, Tuple(species))
    end
    function LennardJones{T}(; ϵ::T, σ::T, rcutoff::T, species::Tuple) where {T}
        new{T}(ϵ, σ, rcutoff, species)
    end
end

get_rcutoff(lj::LennardJones) = lj.rcutoff
get_species(lj::LennardJones) = lj.species

potential_energy(R::AbstractFloat, lj::LennardJones) = 4lj.ϵ * ((lj.σ / R)^12 - (lj.σ / R)^6)
force(R::AbstractFloat, r::SVector{3}, lj::LennardJones) = (24lj.ϵ * (2(lj.σ / R)^12 - (lj.σ / R)^6) / R^2)r
