# InteratomicPotentials Interface

<!--
TODO (Dallas): write explanation of how to create new potentials of various types
example: https://cesmix-mit.github.io/Atomistic.jl/dev/interface/
-->

# Instantiating a Built-In Interatomic Potentials
There are a number of built-in interatomic potentials that are available for users: Born-Mayer-Huggins, Lennard-Jones, Coulomb's Law, and Morse. These potentials are subtypes of `EmpiricalPotential`. There are additional basis potentials, the Spectral Neighbor Analysis Potential and Atomic Cluster Expansion, available in the related package `InteratomicBasisPotentials.jl`. In order to instantiate a potential, you need to provide the necessary potential parameters (see API documentation), radial cutoff, and species (i.e. elements) for which the potential is defined for.

Below is an example for the Lennard-Jones potential:
```julia
ϵ = 0.1034 * u"eV"
σ = 1.0 * u"Å"
rcutoff = 4.0 * σ
species = [:Ar, ]
lj = LennardJones(ϵ, σ, rcutoff, species)
```

In general, empirical potentials handles have the signature:
```julia
MockPotential{T<:AbstractFloat} <: EmpiricalPotential{NamedTuple{(:a, :b, ..., :e)},NamedTuple{(:rcutoff,)}}
# where parameters are :a, :b, ..., :e
```

# Evaluating Interatomic Potentials
There are a number of signatures for evaluating the potential energy, force, and virial stress for the interatomic potentials. For most users it will be most convenient in order to define a configuration in terms of an `AtomsBase.jl` system (i.e., using a `FlexibleSystem`), which requires a list of atoms, boundary conditions, and box information. With such a system, evaluating the potential energy and force are straightforward:

```julia
system = FlexibleSystem(atoms, boundary_conditions, box)
energy, force = energy_and_force(system, lj) # Efficient implementation of energy and force calculation
vs = virial_stress(system, lj)
```

where the abstract signatures are
```julia
energy_and_force(s::AbstractSystem, p::EmpiricalPotential) :: (Unitful.Energy, SVector{n, SVector{3, Unitful.Force}})
virial_stress(s::AbstractSystem, p::EmpiricalPotential) :: SVector{6, Unitful.Energy}
virial(s::AbstractSystem, p::EmpiricalPotential) :: Unitful.Energy
```
where `n` is the number of atoms in the configuration.

Convenience functions exist that return energy and forces alone,
```julia
# Compute both, at the cost of additional neighborlist iterations
potential_energy(s::AbstractSystem, p::EmpiricalPotential) :: Unitful.Energy
force(s::AbstractSystem, p::EmpiricalPotential) :: SVector{n, SVector{3, Unitful.Force}}
```

There are also convenience functions that allow the evaluation of the interatomic potentials as a function of radial distance (where appropriate, for two body potentials), without unitful annotations.
```julia
# Compute both, at the cost of additional neighborlist iterations
potential_energy(R::AbstractFloat, p::EmpiricalPotential) :: AbstractFloat
force(R::AbstractFloat, r::SVector{3}, p::EmpiricalPotential) :: SVector{3}
```

# Manipulating EmpircalPotentials
`MixedPotentials` can be created by algebraically manipulating `EmpiricalPotentials` to create more complex interatomic potentials. To model a mixed configuration of Argon and Xenon potentials, each with their own Lennard-Jones potential and a separate Lennard-Jones potential that defines their interaction, use:
```julia
lj1 = LennardJones(ϵ, σ, rcutoff, [:Ar, ]) # EmpiricalPotential
lj2 = LennardJones(1.35ϵ, 2.0σ, rcutoff, [:Xe, ]) # EmpiricalPotential
lj3 = LennardJones(2ϵ, 1.5σ, rcutoff, [:Ar, :Xe]) # EmpiricalPotential
lj_sum = lj1+lj2+lj3 # MixedPotential
```
When evaluated using an `AtomsBase.jl` system, energies and forces are appropriately calculated.

