"""
    Coulomb{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{()},NamedTuple{(:rcutoff,)}}

The Coulomb, or electrical potential is a simple two-body intermolecular potential that describes the electrical potential energy between two atoms with charge ``q_{1}`` and ``q_{2}``, at a distance ``r``, is given by 

```math
\begin{equation}
V_{LJ}(r; \\epsilon, \\sigma, rcutoff, species) =
\\begin{cases} 
    0 & r > rcutoff \\
    \\frac{k_{e} q_{1} q_{2}}{r} & r < rcutoff.
\\end{cases}
\end{equation}
```
where ``k_{e}`` is known as Coulomb's constant.

Users must supply the electric charges of the atoms, ``q_{1}`` and ``q_{2}``, and the radial cutoff.
"""
struct Coulomb{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{()},NamedTuple{(:rcutoff,)}}
    q₁::T
    q₂::T
    rcutoff::T
    species::Tuple
    function Coulomb(q₁, q₂, rcutoff, species)
        q₁, q₂, rcutoff = promote(austrip(q₁), austrip(q₂), austrip(rcutoff))
        new{typeof(rcutoff)}(q₁, q₂, rcutoff, Tuple(species))
    end
    function Coulomb{T}(; q₁::T, q₂::T, rcutoff::T, species::Tuple) where {T}
        new{T}(q₁, q₂, rcutoff, species)
    end
end

get_rcutoff(c::Coulomb) = c.rcutoff
get_species(c::Coulomb) = c.species

potential_energy(R::AbstractFloat, c::Coulomb) = kₑ * c.q₁ * c.q₂ / R
force(R::AbstractFloat, r::SVector{3}, c::Coulomb) = (kₑ * c.q₁ * c.q₂ / R^3)r
