# [WIP] InteratomicPotentials.jl

[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square")](https://mit-license.org)
[![Build Status](https://github.com/cesmix-mit/InteratomicPotentials.jl/workflows/CI/badge.svg)](https://github.com/cesmix-mit/InteratomicPotentials.jl/actions)
[![codecov](https://codecov.io/gh/cesmix-mit/InteratomicPotentials.jl/branch/main/graph/badge.svg?token=IF6zvl50j9)](https://codecov.io/gh/cesmix-mit/InteratomicPotentials.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://cesmix-mit.github.io/InteratomicPotentials.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesmix-mit.github.io/InteratomicPotentials.jl/dev)

This module implements methods (energies, forces, and virial tensors) for a variety of interatomic potentials, including the SNAP Potential (Thompson et al. 2014, see notes for current limitations).

Project goals:

- Having a defined structure for each potential
  - The potential structure will hold the trainable and nontrainable parameters
  - The potential structure will naturally plug-and-play with fitting and uncertainty quantification codes (PotentialLearning.jl and PotentialUQ.jl)
  - Todo: Develop a struct and interface for Interatomic Potentials that allows the creation of more complex potentials (i.e., sums of potentials or mixed type potentials for systems with multiple elements.)
- Allow for easy use of automatic differentiation through framework.

Right now, this module contains the framework for the following potentials

- Lennard Jones
- SNAP (explicit multi-element support (chem flag) is nominally supported but not currently tested.)

To Dos:

- Improve neighborlist framework for multi-element systems (i.e. allowing potentials to differentiate based on the element of atom i and atom j).
- Finish testing SNAP potential for explicit multi-element flag.
- Create a structure that allows for the complexification and combination of interatomic potentials.
  - Ideally this would allow for the specification of which potential should be used in multi-element systems for each pair of elements.
- Add zbl potential.
- Interface with MDP.
- See issues for more topics.

## Working Example

In order to compute the interatomic energy of a system, or the forces between atoms in a system, the user has to

- 1. define an `AbstractSystem` using ` AtomsBase` and
- 2. construct a potential (subtype of an ArbitraryPotential).

Once these two structures have ben instantiated, the quantity of interest can be computed using the signature `func(system, potential)`.

```julia
# Define an atomic system
atom1     = Atom(element, ( @SVector [1.0, 0.0, 0.0] ) * 1u"Å")
atom2    = Atom(element, ( @SVector [1.0, 0.25, 0.0] ) * 1u"Å")
box = [[1., 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] * 1u"Å"
bcs = [DirichletZero(), Periodic(), Periodic()]
system   = FlexibleSystem(atoms, box , bcs)

ϵ = 1.0
σ = 0.25 * 1u"Å"
rcutoff  = 2.0 * 1u"Å"
lj       = LennardJones(ϵ, σ, rcutoff)           # <: EmpiricalPotential <: ArbitraryPotential
pe       = potential_energy(system, lj)               # <: Float64
f        = force(system, lj)                          # <: Vector{SVector{3, Float64}}
v        = virial(system, lj)                         # <: Float64
v_tensor = virial_stress(system, lj)                  # <: SVector{6, Float64}
```

See "/test/" for further examples.

## Potential Types

All interatomic potentials listed in this project are subtypes of `ArbitraryPotential`. At this point, as of v0.11, there are two branches of potentials: `EmpiricalPotential` and `BasisPotential`. `EmpiricalPotential`s include two-body potentials like `BornMayer`, `LennardJones`. `BasisPotential` require the expansion of the configration of atoms using some basis set of functions and a set of coefficients for each of the members of the expansion in order to calculate energies and forces (i.e. `SNAP`).

```julia
EmpiricalPotential <: ArbitraryPotential
BornMayer <: EmpiricalPotential
LennardJones <: EmpiricalPotential
Coulomb     <: EmpiricalPotential
ZBL         <: EmpiricalPotential

BasisPotential <: ArbitraryPotential
SNAP           <: BasisPotential
```

`BasisPotential` are associated with the functions that evaluate the basis set with parameters necessary for the particular basis expansion used. For example,

```julia
SNAP_parameters = SNAPParams(...)
basis_evaluation = evaluate(system, SNAP_parameters)
gradient_basis_evaulation = evaluate_d(system, SNAP_parameters)
```

For linear basis expansions, the energies and forces are dot products between the potential coefficients and the results of `evaulate` and `evalulate_d`, respectively.
