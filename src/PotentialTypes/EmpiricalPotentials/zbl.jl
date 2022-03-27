## This page is a work in progress. In particular, need to implement switching function that is typically 
## employed within LAMMPS.

struct ZBL <: EmpiricalPotential
    Z_1
    Z_2
    e
    ϵ0
    rcutoff
    species::AbstractVector
end

get_parameters(zbl::ZBL) = Parameter{}(())
set_parameters(p::Parameter{}, zbl::ZBL) = zbl

deserialize_parameters(p::Parameter{()}, zbl::ZBL) = []
serialize_parameters(p::Vector, zbl::ZBL) = Parameter{()}(())

get_hyperparameters(zbl::ZBL) = Parameter{:rcutoff}((zbl.rcutoff,))
set_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = ZBL(zbl.Z_1, zbl.Z_2, zbl.e, zbl.ϵ0, p.rcutoff, zbl.species)

deserialize_hyperparameters(p::Parameter{(:rcutoff,)}, zbl::ZBL) = [p.rcutoff]
serialize_hyperparameters(p::Vector, zbl::ZBL) = Parameter{(:rzblutoff,)}((p[1],))

phi(d::AbstractFloat) = 0.18175 * e^(-3.19980 * d) + 0.50986 * e^(-0.94229 * d) + 0.28022 * e^(-0.40290 * d) + 0.02817 * e^(-0.20162 * d)
dphi_dr(d::AbstractFloat) = -3.19980 * 0.18175 * e^(-3.19980 * d) - 0.94229 * 0.50986 * e^(-0.94229 * d) - 0.40290 * 0.28022 * e^(-0.40290 * d) - 0.20162 * 0.02817 * e^(-0.20162 * d)
############################# Energies ##########################################

function potential_energy(R::AbstractFloat, zbl::ZBL)
    a = 0.8854 * 0.529 / (Z_1^(0.23) + Z_2^(0.23))
    d = R / a
    zbl.Z_1 * zbl.Z_2 * e^2 * phi(d) / (4.0 * π * zbl.ϵ0 * R)
end

############################### Forces ##########################################

function force(R::AbstractFloat, r::SVector{3}, zbl::ZBL)
    a = 0.8854 * 0.529 / (Z_1^(0.23) + Z_2^(0.23))
    d = R / a
    f = zbl.Z_1 * zbl.Z_2 * e^2 * phi(d) / (4.0 * π * zbl.ϵ0 * R^2)
    f += zbl.Z_1 * zbl.Z_2 * e^2 * dphi_dr(d) / (a * 4.0 * π * zbl.ϵ0 * R)
    SVector(f .* r / R)
end

###############################   Gradients  ###################################

# TODO
