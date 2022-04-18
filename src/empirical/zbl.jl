# This implementation is WIP. In particular, need to implement switching function that is typically employed within LAMMPS.
# See https://github.com/cesmix-mit/InteratomicPotentials.jl/issues/11

"""
    ZBL{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{()},NamedTuple{(:rcutoff,)}}

#TODO: DALLAS
The Ziegler-Biersack-Littmark (ZBL) screened nuclear repulsion for describing high-energy collisions between atoms.
"""
struct ZBL{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{()},NamedTuple{(:rcutoff,)}}
    Z₁::Int
    Z₂::Int
    e::T
    rcutoff::T
    species::Tuple
    function ZBL(Z₁, Z₂, e, rcutoff, species)
        Z₁, Z₂, e, rcutoff = promote(Z₁, Z₂, austrip(e), austrip(rcutoff))
        new{typeof(rcutoff)}(Z₁, Z₂, e, rcutoff, Tuple(species))
    end
    function ZBL{T}(; Z₁::Int, Z₂::Int, e::T, rcutoff::T, species::Tuple) where {T}
        new{T}(Z₁, Z₂, e, rcutoff, species)
    end
end

get_rcutoff(zbl::ZBL) = zbl.rcutoff
get_species(zbl::ZBL) = zbl.species

const ϕ_coeffs = [0.18175, 0.50986, 0.28022, 0.02817]
const ϕ_exps = [3.19980, 0.94229, 0.40290, 0.20162]
ϕ(d::AbstractFloat, e::AbstractFloat) = sum(coeff * e^(-exp * d) for (coeff, exp) ∈ zip(ϕ_coeffs, ϕ_exps))
dϕdr(d::AbstractFloat, e::AbstractFloat) = -sum(exp * coeff * e^(-exp * d) for (coeff, exp) ∈ zip(ϕ_coeffs, ϕ_exps))

_zbl_coeff(zbl::ZBL) = kₑ * zbl.Z₁ * zbl.Z₂ * zbl.e^2
_zbl_d(R::AbstractFloat, zbl::ZBL) = R / (0.8854 * 0.529 / (zbl.Z₁^(0.23) + zbl.Z₂^(0.23)))

potential_energy(R::AbstractFloat, zbl::ZBL) = _zbl_coeff(zbl) * ϕ(_zbl_d(R, zbl), zbl.e) / R
force(R::AbstractFloat, r::SVector{3}, zbl::ZBL) = (_zbl_coeff(zbl) * (ϕ(_zbl_d(R, zbl), zbl.e) + _zbl_d(R, zbl) * dϕdr(_zbl_d(R, zbl), zbl.e)) / R^3)r
