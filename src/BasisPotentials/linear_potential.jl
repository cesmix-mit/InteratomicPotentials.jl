"""
    LBasisPotential <: LinearBasisPotential 

Type covering interatomic potentials that produce energies and forces as a linear function of the basis descriptors. Typical examples include vanilla ACE/SNAP. Each struct contains the parameters of the potential, `β', and the basis system being used, `basis'.
"""
struct LBasisPotential{T} <: LinearBasisPotential{NamedTuple{(:β, :β0)}, NamedTuple{()}}
    β::Vector{T}
    β0::Vector{T}
    basis::BasisSystem
end

function LBasisPotential(
    basis :: BasisSystem
)
    return LBasisPotential(zeros(length(basis)), zeros(1), basis)
end

# Get all energies and forces for basis potential

function potential_energy(
    B::Vector{Vector{T}},
    lbp::LBasisPotential{T}
) where T<: Real
    G = compute_global_descriptors(B, lbp)
    return lbp.β0[1] + G ⋅ lbp.β
end

function force(
    dB::Vector{Vector{Vector{T}}},
    lbp::LBasisPotential{T}
) where T<: Real
    f = [[dB_atom_comp' ⋅ lb.β for dB_atom_comp in dB_atom]
         for dB_atom in dB]
    return f
end

function virial_stress(
    W::Matrix{T},
    lbp::LBasisPotential{T}
) where T <: Real 
    return SVector{6}(W * lbp.β)
end

function virial(
    W::Matrix{T},
    lbp::LBasisPotential{T}
) where T<:Real
    return sum(W * lbp.β)
end

# Other useful function

function get_rcutoff(
    lbp::LBasisPotential{T}
) where T<:Real
    return get_rcutoff(lbp.basis)
end

function get_species(
    lbp::LBasisPotential{T}
) where T<:Real
    return get_species(lbp.basis)
end

function get_parameters(
    lbp::LBasisPotential{T}
) where T<:Real
    return (lbp.β0, lbp.β)
end

