"""
    LBasisPotential <: LinearBasisPotential 

Type covering interatomic potentials that produce energies and forces as a linear function of the basis descriptors. Typical examples include vanilla ACE/SNAP. Each struct contains the parameters of the potential, `β', and the basis system being used, `basis'.
"""
struct LBasisPotential{T} <: LinearBasisPotential{NamedTuple{(:β, )}, NamedTuple{()}}
    β :: Vector{T}
    basis :: BasisSystem
end

function LBasisPotential(basis :: BasisSystem)
    LBasisPotential(zeros(length(basis)), basis)
end

get_rcutoff(lbp::LBasisPotential) = get_rcutoff(lbp.basis)
get_species(lbp::LBasisPotential) = get_species(lbp.basis)

potential_energy(B::Vector{Vector{T}}, lbp::LBasisPotential) where T<: Real = sum(dot(lbp.β, b) for b in B)
potential_energy(sys::AtomsBase.AbstractSystem, lbp::LinearBasisPotential) = potential_energy(compute_local_descriptors(sys, lbp.basis), lbp)

force(B::Vector{Matrix{T}}, lbp::LBasisPotential) where T <: Real = SVector{3}.([ b * lbp.β for b in B])
force(s::AtomsBase.AbstractSystem, lbp::LBasisPotential) = force(compute_force_descriptors(s, lbp.basis), lbp)

virial_stress(W::Matrix{T}, lbp::LBasisPotential) where T <: Real = SVector{6}(W * lbp.β)
virial(W::Matrix{T}, lbp::LBasisPotential) where T<:Real = sum(W * lbp.β)

virial_stress(sys::AtomsBase.AbstractSystem, lbp::LBasisPotential) = virial_stress(compute_virial_descriptors(sys, lbp.basis), lbp)
virial(sys::AtomsBase.AbstractSystem, lbp::LBasisPotential) = sum(virial_stress(sys, lbp))

get_parameters(lbp::LBasisPotential) = lbp.β

