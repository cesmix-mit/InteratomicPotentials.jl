"""
    NNBasisPotential <: NeuralNetworkBasisPotential

Definition of the neural network basis potential composed type.
"""
struct NNBasisPotential <: NeuralNetworkBasisPotential{NamedTuple{(:params, )}, NamedTuple{()}}
    nns::Dict
    basis::BasisSystem
end

"""
potential_energy(A::AbstractSystem, p::NNBasisPotential)

`A`: atomic configuration.
`p`: neural network basis potential.

Returns the potential energy of a system using a neural network basis potential.
See 10.1103/PhysRevLett.98.146401, https://fitsnap.github.io/Pytorch.html
"""
function potential_energy(
    A::AbstractSystem,
    nnbp::NNBasisPotential
)
    local_descr = compute_local_descriptors(A, nnbp.basis, T = Float32)
    species = atomic_symbol.(A.particles)
    return sum([nnbp.nns[s](d) for (s, d) in zip(species, local_descr)])[1]
end
