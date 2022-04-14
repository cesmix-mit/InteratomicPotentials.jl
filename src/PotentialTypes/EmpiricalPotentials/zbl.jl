## This page is a work in progress. In particular, need to implement switching function that is typically 
## employed within LAMMPS.

############################## ZBL ###################################
struct ZBL <: EmpiricalPotential
    Z_1
    Z_2
    e
    rcutoff
    species::AbstractVector
end
function ZBL(Z_1::Integer, Z_2::Integer, e::Unitful.Charge, rcutoff::Unitful.Length, species::AbstractVector)
    ZBL(Z_1, Z_2, austrip(e), austrip(rcutoff), species)
end

get_parameters(zbl::ZBL) = Parameter{}(())
set_parameters(p::Parameter{}, zbl::ZBL) = zbl

deserialize_parameters(p::Parameter{()}, zbl::ZBL) = []
serialize_parameters(p::Vector, zbl::ZBL) = Parameter{()}(())

get_hyperparameters(zbl::ZBL) = Parameter{:rcutoff}((zbl.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = ZBL(zbl.Z_1, zbl.Z_2, zbl.e, p.rcutoff, zbl.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = [p.rcutoff]
serialize_hyperparameters(p::Vector, zbl::ZBL) = Parameter{(:rzblutoff,)}((p[1],))

_ϕ(d::AbstractFloat, e::AbstractFloat) = 0.18175 * e^(-3.19980 * d) + 0.50986 * e^(-0.94229 * d) + 0.28022 * e^(-0.40290 * d) + 0.02817 * e^(-0.20162 * d)
_dϕdr(d::AbstractFloat, e::AbstractFloat) = -3.19980 * 0.18175 * e^(-3.19980 * d) - 0.94229 * 0.50986 * e^(-0.94229 * d) - 0.40290 * 0.28022 * e^(-0.40290 * d) - 0.20162 * 0.02817 * e^(-0.20162 * d)
_zbl_coeff(zbl::ZBL) = kₑ * zbl.Z_1 * zbl.Z_2 * zbl.e^2
_zbl_d(R::AbstractFloat, zbl::ZBL) = R / (0.8854 * 0.529 / (zbl.Z_1^(0.23) + zbl.Z_2^(0.23)))

potential_energy(R::AbstractFloat, zbl::ZBL) = _zbl_coeff(zbl) * _ϕ(_zbl_d(R, zbl), zbl.e) / R
force(R::AbstractFloat, r::SVector{3}, zbl::ZBL) = (_zbl_coeff(zbl) * (_ϕ(_zbl_d(R, zbl), zbl.e) + _zbl_d(R, zbl) * _dϕdr(_zbl_d(R, zbl), zbl.e)) / R^3)r
