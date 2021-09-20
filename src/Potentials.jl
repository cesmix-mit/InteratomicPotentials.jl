################################################################################
#
#    Module Potentials.jl
#
#    This module will implement a variety of interatomic potentials and defines abstract types
#    that allow these potentials to be used in other packages (PotentialLearning.jl, PotentialUQ.jl, 
#    Atomistic.jl, etc...). 
#       'implementation' means:
#           1. Having a defined structure for each potential
#               1.1 The structure holds all of the necessary parameters for evaluating energies, forces, ...
#               1.2 The potential structure should expose the trainable and nontrainable parameters (necessary for learning).
#           2. Having a method to get the potential energy of given configuration, as
#               defined by that potential.
#           3. Having a method to produce the force of a given configuration, as defined
#               by that potential.
#           4. Having a method to produce the stresses of a given configuration, as defined
#               by that potential.
#           5. (For inference) Having a method to produce the gradient of each of the above methods with 
#               respect to the potential parameters.
#
#       Right now, this module contains the framework for the following potentials
#           1. Lennard Jones
#           2. Born - Mayer 
#           3. Coulomb
#           4. ZBL
#           5. GaN (special mixed type, more of an application than a true potential)
#           6. SNAP 

#
#       To do:
#           1. Improve configuration implementation
#           2. Implement ACE potentials (with ACE.jl)
#
#
################################################################################
module Potentials

using Base: Float64
using LAMMPS
using LinearAlgebra
using Plots

include("Utilities/utils.jl")
include("Configurations/config.jl")
include("IO/io.jl")
include("PotentialTypes/types.jl")
include("MD/md.jl")






end