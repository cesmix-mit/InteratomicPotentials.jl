export LinearBasisPotential, NNBasisPotential
include("linear_potential.jl")
include("nn_potential.jl")

################################################################################
# InteratomicPotentials API implmentations for basis potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::LinearBasisPotential)
    B = get_local_descriptors(A, p)
    dB = get_force_descriptors(A, p)
    e = sum(B) ⋅ p.coefficients # * ENERGY_UNIT
    f = [SVector{3}(di * p.coefficients) for di in dB] # * FORCE_UNIT
    (; e, f)
end

function energy_and_force(A::AbstractSystem, p::NNBasisPotential)
    B = sum(get_local_descriptors(A, p))
    dB = get_force_descriptors(A, p)

    e = p.nn(B)
    e, dnndb = withgradient(p.nn, b)
    dbdr = get_force_descriptors(A, p.basis)
    f = [SVector{3}.([-dnndb ⋅ dbdr[atom][coor,:] for coor in 1:3])
         for atom in 1:length(dbdr)]
    (; e, f)
end

# function virial_stress(A::AbstractSystem, p::BasisPotential)
#     sum(SVector{6}(di ⋅ p.coefficients) for di in evaluate_basis_v(A, p.basis_params)) * ENERGY_UNIT
# end