# [WIP] InteratomicPotentials.jl
This module implements methods (energies, forces, and virial tensors) for a variety of interatomic potentials, including the SNAP Potential (Thompson et al. 2014). 

Project goals:
- Having a defined structure for each potential
    - The potential structure will hold the trainable and nontrainable parameters 
    - The potential structure will naturally plug-and-play with fitting and uncertainty quantification codes (PotentialLearning.jl and PotentialUQ.jl)
- Have a defined structure for each configuration of atoms.
    - A configuration of atom holds identifying information for a group of atoms (atom masses, names, positions, velocities, etc..).
    - Energy, force, and virial methods are implemented for these configurations.
    - Configurations have enough information to allow for use with LAMMPS.jl (and potentially other codes in the future). 
    - Input and output of atomic data is supported through configurations. 
- Allow for easy use of automatic differentiation through framework.
- Provide fledgling support for molecular dynamics codes (potentially LAMMPS.jl, NBodySimulator.jl, and others...)

Right now, this module contains the framework for the following potentials
- Lennard - Jones
- Born - Mayer 
- Coulomb
- GaN (special mixed type, more of an application than a true potential)
- SNAP 

## Working Example
Load a configuration of Argon atoms in a periodic system with Lennard-Jones potential and units ([LAMMPS DATA File](https://docs.lammps.org/2001/data_format.html)).
```julia
r        = load_lammps_data(file_path; atom_names = [:Ar],              
                                radii = [1.0], weights = [1.0], 
                                boundaries = ["p", "p", "p"], units = "lj")  
                                # typeof(r) <: Potentials.Configuration
ϵ        = 1
σ        = 1
lj       = LennardJones(ϵ, σ)                    # <: EmpiricalPotential <: ArbitraryPotential
pe       = potential_energy(r, lj)               # <: Real                    
f        = Potentials.force(r, lj)               # <: Vector{Real}(, 3)
v        = Potentials.virial(r, lj)              # <: Real
v_tensor = Potentials.virial_stress(r, lj)       # <: Vector{Real}(, 6)
```
See "/test/" for further examples.


To do:
- Make usable with PotentialLearning.jl and PotentialUQ.jl
- Finish implemented gradients of energies, forces, stresses w.r.t. potential parameters.
- Improve input and output support.
- Improve molecular dynamics support.
