# [WIP] InteratomicPotentials.jl
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

## Working Example
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
pe       = potential_energy(r, lj)               # <: Float64                   
f        = Potentials.force(r, lj)               # <: SVector{2, SVector{3, Float64}}
v        = Potentials.virial(r, lj)              # <: Float64
v_tensor = Potentials.virial_stress(r, lj)       # <: SVector{6, Float64}
```
See "/test/" for further examples.


