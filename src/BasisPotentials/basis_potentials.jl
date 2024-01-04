################################################################################
# InteratomicPotentials API implementations for basis potentials
################################################################################

# Common basis potential functions for computing energy, forces, and virials

function compute_global_descriptors(
    B::Vector{Vector{T}},
    bp::BasisPotential
) where T<: Real
    return sum(B)
end

function potential_energy(
    s::AbstractSystem,
    bp::BasisPotential
)
    B = compute_local_descriptors(s, bp.basis)
    e = potential_energy(B, bp)
    return e
end

function force(
    s::AbstractSystem,
    bp::BasisPotential
)
    dB = compute_force_descriptors(s, bp.basis)
    f = force(dB, bp) # * FORCE_UNIT
    return f
end

function energy_and_force(
    s::AbstractSystem,
    bp::BasisPotential
)
    B = compute_local_descriptors(s, bp.basis)
    e = potential_energy(B, bp) # * ENERGY_UNIT
    dB = compute_force_descriptors(s, bp.basis)
    f = force(dB, bp) # * FORCE_UNIT
    (; e, f)
end

function virial_stress(
    s::AbstractSystem,
    bp::BasisPotential
)
    W = compute_virial_descriptors(s, bp)
    v = virial_stress(W, bp)
    return v
end

function virial(
    s::AbstractSystem,
    lbp::BasisPotential
)
    return sum(virial_stress(s, lbp))
end

# Specific basis potential functions for computing energy, forces, and virials

include("linear_potential.jl")
include("nn_potential.jl")
