"""
    NNBasisPotential <: NeuralNetworkBasisPotential

Definition of the neural network basis potential composed type.
"""
struct NNBasisPotential <: NeuralNetworkBasisPotential{NamedTuple{(:params, )}, NamedTuple{()}}
    nns :: Dict
    basis :: BasisSystem
end

"""
potential_energy(A::AbstractSystem, p::NNBasisPotential)

`c`: atomic configuration.
`p`: neural network basis potential.

Returns the potential energy of a system using a neural network basis potential.
See 10.1103/PhysRevLett.98.146401, https://fitsnap.github.io/Pytorch.html
"""
function potential_energy(
    c::Configuration,
    nnbp::NNBasisPotential
)
    local_descr = get_values(get_local_descriptors(c))
    species = atomic_symbol.(get_system(c).particles)
    return sum([nnbp.nns[s](d) for (s, d) in zip(species, local_descr)])[1]
end
