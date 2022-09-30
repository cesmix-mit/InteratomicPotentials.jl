"""
    LinearBasisPotential <: BasisPotential 

Type covering interatomic potentials that produce energies and forces as a linear function of the basis descriptors. Typical examples include vanilla ACE/SNAP. Each struct contains the parameters of the potential, `β', and the basis system being used, `basis'.
"""
struct LinearBasisPotential{T} <: BasisPotential{NamedTuple{(:β, )}, NamedTuple{()}}
    β :: Vector{T}
    basis :: BasisSystem
end

function LinearBasisPotential(basis :: BasisSystem)
    LinearBasisPotential(zeros(length(basis)), basis)
end

get_rcutoff(lbp::LinearBasisPotential) = get_rcutoff(lbp.basis)
get_species(lbp::LinearBasisPotential) = get_species(lbp.basis)

potential_energy(B::Vector{Vector{T}}, lbp::LinearBasisPotential) where T<: Real = sum(dot(lbp.β, b) for b in B)
potential_energy(sys::AtomsBase.AbstractSystem, lbp::LinearBasisPotential) = potential_energy(get_local_descriptors(sys, lbp.basis), lbp)

force(B::Vector{Matrix{T}}, lbp::LinearBasisPotential) where T <: Real = SVector{3}.([ b * lbp.β for b in B])
force(s::AtomsBase.AbstractSystem, lbp::LinearBasisPotential) = force(get_force_descriptors(s, lbp.basis), lbp)

virial_stress(W::Matrix{T}, lbp::LinearBasisPotential) where T <: Real = SVector{6}(W * lbp.β)
virial(W::Matrix{T}, lbp::LinearBasisPotential) where T<:Real = sum(W * lbp.β)

virial_stress(sys::AtomsBase.AbstractSystem, lbp::LinearBasisPotential) = virial_stress(get_virial_descriptors(sys, lbp.basis), lbp)
virial(sys::AtomsBase.AbstractSystem, lbp::LinearBasisPotential) = sum(virial_stress(sys, lbp))

get_parameters(lbp::LinearBasisPotential) = lbp.β

