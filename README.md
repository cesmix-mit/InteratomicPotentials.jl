# [WIP] Potentials.jl
This module will implement a variety of interatomic potentials. 

Here, 'implementation' means:
- Having a defined structure for each potential
    - The potential structure will hold the trainable and nontrainable parameters 
- Having a method to get the potential energy of a given configuration, as defined by that potential.
- Having a method to produce the force of a given configuration, as defined by that potential.
- Having a method to produce the stresses of a given configuration, as defined by that potential.

Right now, this module contains the framework for the following potentials
- Lennard - Jones
- Born - Mayer 
- Coulomb
- GaN (special mixed type, more of an application than a true potential)
- SNAP 

To do:
- Make usable with PotentialLearning.jl and PotentialUQ.jl
- Finish implemented gradients of energies, forces, stresses w.r.t. potential parameters.
