## This page is a work in progress. In particular, need to implement switching function that is typically 
## employed within LAMMPS.

############################## ZBL ###################################
struct ZBL{T<:AbstractFloat} <: EmpiricalPotential
    Z₁::Int
    Z₂::Int
    e::T
    rcutoff::T
    species::Vector{Symbol}
end
function ZBL(Z₁::Integer, Z₂::Integer, e::Unitful.Charge, rcutoff::Unitful.Length, species::AbstractVector{Symbol})
    ZBL(Z₁, Z₂, austrip(e), austrip(rcutoff), collect(species))
end

get_parameters(zbl::ZBL) = Parameter{}(())
set_parameters(p::Parameter{}, zbl::ZBL) = zbl

deserialize_parameters(p::Parameter{()}, zbl::ZBL) = []
serialize_parameters(p::Vector, zbl::ZBL) = Parameter{()}(())

get_hyperparameters(zbl::ZBL) = Parameter{:rcutoff}((zbl.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = ZBL(zbl.Z₁, zbl.Z₂, zbl.e, p.rcutoff, zbl.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = [p.rcutoff]
serialize_hyperparameters(p::Vector, zbl::ZBL) = Parameter{(:rzblutoff,)}((p[1],))

const _ϕ_coeffs = [0.18175, 0.50986, 0.28022, 0.02817]
const _ϕ_exps = [3.19980, 0.94229, 0.40290, 0.20162]
_ϕ(d::AbstractFloat, e::AbstractFloat) = sum(coeff * e^(-exp * d) for (coeff, exp) ∈ zip(_ϕ_coeffs, _ϕ_exps))
_dϕdr(d::AbstractFloat, e::AbstractFloat) = -sum(exp * coeff * e^(-exp * d) for (coeff, exp) ∈ zip(_ϕ_coeffs, _ϕ_exps))

_zbl_coeff(zbl::ZBL) = kₑ * zbl.Z₁ * zbl.Z₂ * zbl.e^2
_zbl_d(R::AbstractFloat, zbl::ZBL) = R / (0.8854 * 0.529 / (zbl.Z₁^(0.23) + zbl.Z₂^(0.23)))

potential_energy(R::AbstractFloat, zbl::ZBL) = _zbl_coeff(zbl) * _ϕ(_zbl_d(R, zbl), zbl.e) / R
force(R::AbstractFloat, r::SVector{3}, zbl::ZBL) = (_zbl_coeff(zbl) * (_ϕ(_zbl_d(R, zbl), zbl.e) + _zbl_d(R, zbl) * _dϕdr(_zbl_d(R, zbl), zbl.e)) / R^3)r
