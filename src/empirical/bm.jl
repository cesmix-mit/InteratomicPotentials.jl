##############################   Born-Mayer  ###################################
struct BornMayer{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:A, :ρ)},NamedTuple{(:rcutoff,)}}
    A::T
    ρ::T
    rcutoff::T
    species::Tuple
    function BornMayer(A, ρ, rcutoff, species)
        A, ρ, rcutoff = promote(austrip(A), austrip(ρ), austrip(rcutoff))
        new{typeof(rcutoff)}(A, ρ, rcutoff, Tuple(species))
    end
    function BornMayer{T}(; A::T, ρ::T, rcutoff::T, species::Tuple) where {T}
        new{T}(A, ρ, rcutoff, species)
    end
end

get_rcutoff(bm::BornMayer) = bm.rcutoff
get_species(bm::BornMayer) = bm.species

potential_energy(R::AbstractFloat, bm::BornMayer) = bm.A * exp(-R / bm.ρ)
force(R::AbstractFloat, r::SVector{3}, bm::BornMayer) = (bm.A * exp(-R / bm.ρ) / (bm.ρ * R))r
