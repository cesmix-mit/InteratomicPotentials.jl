## This page is a work in progress. In particular, need to implement switching function that is typically 
## employed within LAMMPS.

############################## ZBL ###################################
struct ZBL{T<:AbstractFloat} <: EmpiricalPotential
    Z₁::Int
    Z₂::Int
    e::T
    rcutoff::T
    species::Tuple
    function ZBL(Z₁, Z₂, e, rcutoff, species)
        Z₁, Z₂, e, rcutoff = promote(Z₁, Z₂, austrip(e), austrip(rcutoff))
        new{typeof(rcutoff)}(Z₁, Z₂, e, rcutoff, Tuple(species))
    end
end

get_rcutoff(zbl::ZBL) = zbl.rcutoff
get_species(zbl::ZBL) = zbl.species

get_parameters(::ZBL) = Parameter{}()
set_parameters(zbl::ZBL, ::Parameter{}) = zbl

serialize_parameters(zbl::ZBL) = collect(get_parameters(zbl))
deserialize_parameters(zbl::ZBL, p::AbstractVector) = set_parameters(zbl, Parameter{}(p))

get_hyperparameters(zbl::ZBL) = Parameter{(:rcutoff,)}((zbl.rcutoff,))
set_hyperparameters(zbl::ZBL, p::Parameter{(:rcutoff,)}) = ZBL(zbl.Z₁, zbl.Z₂, zbl.e, p.rcutoff, zbl.species)

serialize_hyperparameters(zbl::ZBL) = collect(get_hyperparameters(zbl))
deserialize_hyperparameters(zbl::ZBL, p::AbstractVector) = set_hyperparameters(zbl, Parameter{(:rcutoff,)}(p))

const _ϕ_coeffs = [0.18175, 0.50986, 0.28022, 0.02817]
const _ϕ_exps = [3.19980, 0.94229, 0.40290, 0.20162]
_ϕ(d::AbstractFloat, e::AbstractFloat) = sum(coeff * e^(-exp * d) for (coeff, exp) ∈ zip(_ϕ_coeffs, _ϕ_exps))
_dϕdr(d::AbstractFloat, e::AbstractFloat) = -sum(exp * coeff * e^(-exp * d) for (coeff, exp) ∈ zip(_ϕ_coeffs, _ϕ_exps))

_zbl_coeff(zbl::ZBL) = kₑ * zbl.Z₁ * zbl.Z₂ * zbl.e^2
_zbl_d(R::AbstractFloat, zbl::ZBL) = R / (0.8854 * 0.529 / (zbl.Z₁^(0.23) + zbl.Z₂^(0.23)))

potential_energy(R::AbstractFloat, zbl::ZBL) = _zbl_coeff(zbl) * _ϕ(_zbl_d(R, zbl), zbl.e) / R
force(R::AbstractFloat, r::SVector{3}, zbl::ZBL) = (_zbl_coeff(zbl) * (_ϕ(_zbl_d(R, zbl), zbl.e) + _zbl_d(R, zbl) * _dϕdr(_zbl_d(R, zbl), zbl.e)) / R^3)r
