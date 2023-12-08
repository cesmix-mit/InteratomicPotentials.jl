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
    bp::BasisPotential;
    T = Float64
)
    B = compute_local_descriptors(s, bp.basis, T = T)
    e = potential_energy(B, bp)
    return e
end

function force(
    s::AbstractSystem,
    bp::BasisPotential;
    T = Float64
)
    dB = compute_force_descriptors(s, bp.basis, T = T)
    f = force(dB, bp) # * FORCE_UNIT
    return f
end

function energy_and_force(
    s::AbstractSystem,
    bp::BasisPotential;
    T = Float64
)
    B = compute_local_descriptors(s, bp.basis, T = T)
    e = potential_energy(B, bp) # * ENERGY_UNIT
    dB = compute_force_descriptors(s, bp.basis, T = T)
    f = force(dB, bp) # * FORCE_UNIT
    (; e, f)
end

function virial_stress(
    s::AbstractSystem,
    bp::BasisPotential;
    T = Float64
)
    W = compute_virial_descriptors(s, bp, T = T)
    v = virial_stress(W, bp)
    return v
end

function virial(
    s::AbstractSystem,
    lbp::BasisPotential;
    T = Float64
)
    return sum(virial_stress(s, lbp, T = T))
end

# Specific basis potential functions for computing energy, forces, and virials

include("linear_potential.jl")
include("nn_potential.jl")
