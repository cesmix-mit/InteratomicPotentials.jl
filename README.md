# InteratomicPotentials.jl

<!--
TODO (Dallas): The description here is very out of date.
-->

[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square")](https://mit-license.org)
[![Build Status](https://github.com/cesmix-mit/InteratomicPotentials.jl/workflows/CI/badge.svg)](https://github.com/cesmix-mit/InteratomicPotentials.jl/actions)
[![codecov](https://codecov.io/gh/cesmix-mit/InteratomicPotentials.jl/branch/main/graph/badge.svg?token=IF6zvl50j9)](https://codecov.io/gh/cesmix-mit/InteratomicPotentials.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesmix-mit.github.io/InteratomicPotentials.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesmix-mit.github.io/InteratomicPotentials.jl/dev)

This repository implements some basic language and syntax for manipulating interatomic potentials in Julia. The primary purpose of this package is to design a flexible package for use with data-driven and parameter-fitted interatomic potentials. This package is also being designed in order to allow users to define custom potentials and forces for use in molecular dynamics.

This package is part of the CESMIX molecular modeling suite. This package is also intended to be used with Atomistic.jl (for molecular dynamics, with Molly.jl), InteratomicBasisPotentials.jl (for machine learning potentials like SNAP and ACE), and  PotentialLearning.jl (for fitting potentials from data).

This package is a work in progress. 

## Working Example

In order to compute the interatomic energy of a system, or the forces between atoms in a system, the user has to

- 1. define an `AbstractSystem` using ` AtomsBase` and
- 2. construct a potential (subtype of an ArbitraryPotential).

Once these two structures have ben instantiated, the quantity of interest can be computed using the signature `func(system, potential)`.

First, let's create a configuration:
```julia
using AtomsBase, Unitful, UnitfulAtomic
# Define an atomic system
element = :Ar
atom1     = Atom(element, ( @SVector [1.0, 0.0, 0.0] ) * 1u"Å")
atom2    = Atom(element, ( @SVector [1.0, 0.25, 0.0] ) * 1u"Å")
box = [[1., 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] * 1u"Å"
bcs = [DirichletZero(), Periodic(), Periodic()]
system   = FlexibleSystem(atoms, box , bcs)
```
Now we can define the parameters of our interatomic potential:
```julia
ϵ = 1.0 * 1u"eV"
σ = 0.25 * 1u"Å"
rcutoff  = 2.25 * 1u"Å"
lj       = LennardJones(ϵ, σ, rcutoff, [element])           # <: EmpiricalPotential <: AbstractPotential
```

Now we can compute a variety of quantities of the system:
```julia
pe       = potential_energy(system, lj)               # <: Float64 (Hartree)
f        = force(system, lj)                          # <: Vector{SVector{3, Float64}} (Hartree/Bohr)
v        = virial(system, lj)                         # <: Float64 (Hartree)
v_tensor = virial_stress(system, lj)                  # <: SVector{6, Float64} (Hartree)
```

When computing the force, the energy is already available. A convenience implementation that returns both quantities is given by:
```julia 
pe, f    = energy_and_force(system, lj)
```

See "/test/" for further examples.

## Utility functions
There are a growing number of features designed to allow handle of potential parameter easier. For example, one can retrieve the parameters of a potential via:
```julia
get_rcutoff(lj) # Gets radial cutoff (here: 2.25 * 1u"Angstrom")
get_species(lj) # Returns the species the potential is defined for (here: [:Ar])
get_parameters(lj) # Returns the parameters (here: [ϵ, σ])
set_parameters(lj, (ϵ = 2.0 * 1u"eV", σ = 1.0 * 1u"Å")) # Set parameters (returns a new potential)
```

## Potential Types

All interatomic potentials listed in this project are subtypes of `ArbitraryPotential`. At this point, as of v2.0, there are two branches of potentials: `EmpiricalPotential` and `MixedPotentials`. A sister package, ```julia InteratomicBasisPotentials.jl``` defines a potential called `BasisPotential`, see that package for additional details.

`EmpiricalPotential`s include two-body potentials like `BornMayer`, `LennardJones`. `MixedPotential` is a convenience type for allowing the linear combination of potentials. An example would be:
```julia
lj1 = LennardJones(1.0 * u"eV", 1.0 * u"Angstrom", 2.5 * u"Angstrom", [:Ar]) # Ar-Ar Interactions
lj2 = LennardJones(1.0 * u"eV", 1.5 * u"Angstrom", 3.0 * u"Angstrom", [:Xe]) # Xe-Xe Interactions
lj3 = LennardJones(1.5 * u"eV", 1.3 * u"Angstrom", 2.5 * u"Angstrom", [:Ar, :Xe]) # Ar-Xe Interactions
lj = lj1 + lj2 + lj3# Potential defined for all interactions in an Ar-Xe system.
```

```julia
EmpiricalPotential <: AbstractPotential
BornMayer <: EmpiricalPotential
LennardJones <: EmpiricalPotential
Coulomb     <: EmpiricalPotential
ZBL         <: EmpiricalPotential

LinearCombinationPotential <: MixedPotential

# See InteratomicBasisPotentials.jl
BasisPotential <: AbstractPotential
SNAP           <: BasisPotential
ACE            <: BasisPotential
```

