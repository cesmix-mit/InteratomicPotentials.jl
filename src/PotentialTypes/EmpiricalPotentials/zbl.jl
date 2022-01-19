## This page is a work in progress. In particular, need to implement switching function that is typically 
## employed within LAMMPS.

struct ZBL <: EmpiricalPotential
    Z_1 :: AbstractFloat
    Z_2 :: AbstractFloat
    e   :: AbstractFloat
    ϵ0  :: AbstractFloat
    rcutoff :: AbstractFloat
end

get_trainable_params(zbl::ZBL) = Parameter{}(())
get_nontrainable_params(zbl::ZBL) = Parameter{:Z_1, :Z_2, :e, :ϵ0, :rcutoff}((Z_1, Z_2, e, ϵ0, rcutoff))

phi(d::AbstractFloat) =  0.18175 * e^(-3.19980 * d) + 0.50986 * e^(-0.94229 * d) + 0.28022 * e^(-0.40290 * d) + 0.02817 * e^(-0.20162 * d)
dphi_dr(d::AbstractFloat) = -3.19980 * 0.18175 * e^(-3.19980 * d) - 0.94229 * 0.50986 * e^(-0.94229 * d) - 0.40290 * 0.28022 * e^(-0.40290 * d) - 0.20162 * 0.02817 * e^(-0.20162 * d)
############################# Energies ##########################################

function potential_energy(R::AbstractFloat, zbl::ZBL)
    a = 0.8854 * 0.529 / (Z_1^(0.23) + Z_2^(0.23))
    d = R / a
    zbl.Z_1 * zbl.Z_2 * e^2 * phi(d) / (4.0 * π * zbl.ϵ0 * R)
end

############################### Forces ##########################################

force(r::SVector{3,<:AbstractFloat}, zbl::ZBL) = force(r, norm(r), zbl)

function force(R::AbstractFloat, r::SVector{3,<:AbstractFloat}, zbl::ZBL)
    a = 0.8854 * 0.529 / (Z_1^(0.23) + Z_2^(0.23))
    d = R / a
    f = zbl.Z_1 * zbl.Z_2 * e^2 * phi(d) / (4.0 * π * zbl.ϵ0 * R^2) 
    f += zbl.Z_1 * zbl.Z_2 * e^2 * dphi_dr(d) / (a * 4.0 * π * zbl.ϵ0 * R) 
    SVector(f .* r / R)
end

# ##############################   Gradients  ###################################
# TO DO