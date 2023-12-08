"""
    LBasisPotential <: LinearBasisPotential 

Type covering interatomic potentials that produce energies and forces as a linear function of the basis descriptors. Typical examples include vanilla ACE/SNAP. Each struct contains the parameters of the potential, `β', and the basis system being used, `basis'.
"""
struct LBasisPotential{T} <: LinearBasisPotential{NamedTuple{(:β, :β0)}, NamedTuple{()}}
    β
    β0
    basis
end

function LBasisPotential(basis :: BasisSystem; T = Float64)
    return LBasisPotential{T}(zeros(T, length(basis)), zeros(T, 1), basis)
end

# Get all energies and forces for basis potential

function potential_energy(
    B::Vector{Vector{T}},
    lbp::LBasisPotential
) where T<: Real
    G = compute_global_descriptors(B, lbp)
    return lbp.β0[1] + dot(G, lbp.β)
end

function force(
    dB::Vector{Matrix{T}},
    lbp::LBasisPotential
) where T<: Real

#    force_descriptors = [reduce(vcat, get_values(get_force_descriptors(dsi)) ) for dsi in ds]
#    return vcat([lb.β0[1] .+  dB' * lb.β
#                 for dB in [reduce(hcat, fi)
#                 for fi in force_descriptors]]...)

#    G = compute_global_descriptors(B, lbp)
#    return vcat([lb.β0[1] .+  dB' * lb.β for dB in [reduce(hcat, fi) for fi in force_descriptors]]...)
#    
    #force(B::Vector{Matrix{T}}, lbp::LBasisPotential) where T <: Real = SVector{3}.([ b * lbp.β for b in B])
end

function virial_stress(
    W::Matrix{T},
    lbp::LBasisPotential
) where T <: Real 
    return SVector{6}(W * lbp.β)
end

function virial(
    W::Matrix{T},
    lbp::LBasisPotential
) where T<:Real
    return sum(W * lbp.β)
end

# Other useful function

function get_rcutoff(
    lbp::LBasisPotential
)
    return get_rcutoff(lbp.basis)
end

function get_species(
    lbp::LBasisPotential
)
    return get_species(lbp.basis)
end

function get_parameters(
    lbp::LBasisPotential
)
    return (lbp.β0, lbp.β)
end

