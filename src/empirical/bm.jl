"""
    BornMayer{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:A, :ρ)},NamedTuple{(:rcutoff,)}}

The Born-Mayer-Huggins potential is a two-body intermolecular potential with five typical parameters, primarily describing the interaction of neutral atoms. Formally, the interaction between two atoms at a distance, ``r``, is given by 

```math
\begin{aligned}
V_{LJ}(r; \epsilon, \sigma, rcutoff, species) =
    \begin{cases} 
    0 & r > rcutoff \\
    A e^{\frac{\sigma-r}{\rho}} + \left( \frac{C}{r} \right)^8 - \left( \frac{D}{r} \right) & r < rcutoff.
    \end{cases}
\end{aligned}
```

Users must supply five parameters, ``A`` (units of energy), ``\\sigma`` (units of distance), ``\\rho`` (units of distance), ``C`` ``\\&`` ``D`` (units of energy ``\\times`` units of distance), and radial cutoff (units of distance).
"""
struct BornMayer{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:A, :ρ, :σ, :C, :D)},NamedTuple{(:rcutoff,)}}
    A::T
    ρ::T
    σ::T
    C::T
    D::T 
    rcutoff::T
    species::Tuple
    function BornMayer(A, ρ, σ, C, D, rcutoff, species)
        A, ρ, σ, C, D, rcutoff = promote(austrip(A), austrip(ρ), austrip(σ), austrip(C), austrip(D), austrip(rcutoff))
        new{typeof(rcutoff)}(A, ρ, σ, C, D, rcutoff, Tuple(species))
    end
    function BornMayer{T}(; A::T, ρ::T, σ::T, C::T, D::T, rcutoff::T, species::Tuple) where {T}
        new{T}(A, ρ, σ, C, D, rcutoff, species)
    end
end

get_rcutoff(bm::BornMayer) = bm.rcutoff
get_species(bm::BornMayer) = bm.species

potential_energy(R::AbstractFloat, bm::BornMayer) = bm.A * exp((bm.σ - R) / bm.ρ) + (bm.C / R)^6 - (bm.D/R)^8
force(R::AbstractFloat, r::SVector{3}, bm::BornMayer) = (bm.A * exp((bm.σ-R) / bm.ρ)/ (bm.ρ*R) + ((bm.C/R)^6 - (bm.D/R)^8)/R^2)r
