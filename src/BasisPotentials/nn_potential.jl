"""
    NNBasisPotential <: NeuralNetworkBasisPotential

Definition of the neural network basis potential composed type.
"""
struct NNBasisPotential <: NeuralNetworkBasisPotential{NamedTuple{(:params, )}, NamedTuple{()}}
    nns::Dict
    basis::BasisSystem
end

function potential_energy(
    B::Vector{Vector{T}},
    nnbp::NNBasisPotential
) where T<: Real
    species = keys(nnbp.nns) # atomic_symbol.(A.particles)
    return sum([nnbp.nns[s](d) for (s, d) in zip(species, B)])[1]
end

function force(
    dB::Vector{Matrix{T}},
    nnbp::NNBasisPotential
) where T<: Real
    # TODO
    return
end

function virial_stress(
    W::Matrix{T},
    nnbp::NNBasisPotential
) where T <: Real 
    # TODO
    return
end

function virial(
    W::Matrix{T},
    nnbp::NNBasisPotential
) where T<:Real
    # TODO
    return
end
