# The InteratomicBasisPotentials.jl API specification.
import InteratomicPotentials: energy_and_force, virial_stress, ENERGY_UNIT, FORCE_UNIT
import InteratomicPotentials: get_rcutoff, get_species 
import InteratomicPotentials: potential_energy, force, virial
import InteratomicPotentials: compute_local_descriptors, compute_force_descriptors, compute_virial_descriptors, compute_all_descriptors
export energy_and_force, potential_energy, force, virial, virial_stress
export get_rcutoff, get_species, get_rpi
export compute_local_descriptors, compute_force_descriptors, compute_virial_descriptors, compute_all_descriptors

import Base.length
export length
"""
    energy_and_force(s::AbstractSystem, p::AbstractPotential)::NamedTuple{(:e, :f), Tuple{Unitful.Energy,Vector{SVector{3, Unitful.Force}}}}

Calculate the unit-annotated potential energy of a system and the force acting on each particle in a system using the provided interatomic potential.
This combined function is offered because it is usually more efficient to calculate both properties simultaneously.
"""
function energy_and_force end
"""
    potential_energy(s::AbstractSystem, p::AbstractPotential)::Unitful.Energy

Calculate the unit-annotated potential energy of a system using the provided interatomic potential.
The default implementation uses the `:e` property of `energy_and_force(s,p)`.
"""
function potential_energy end
"""
    force(s::AbstractSystem, p::AbstractPotential)::Vector{SVector{3, Unitful.Force}}

Calculate the unit-annotated force acting on each particle in a system using the provided interatomic potential.
The default implementation uses the `:f` property of `energy_and_force(s,p)`.
"""
function force end
"""
    virial_stress(s::AbstractSystem, p::AbstractPotential)::SVector{6,Unitful.Energy}

Calculate the unit-annotated virial stress tensor of a system, officially calculated as the sum of radial-force outerproducts: `` \\sum r_{ij} \\bigotimes F_{ij}``, only returns the unique lower-diagonal components.
"""
function virial_stress end
"""
    virial(s::AbstractSystem, p::AbstractPotential)::Unitful.Energy

Calculate the unit-annotated virial of a system, officially calculated as the trace contraction of the sum of radial-force outerproducts: ``tr\\left( \\sum r_{ij} \\bigotimes F_{ij} \\right)``
"""
function virial end
"""
    get_rcutoff(p::AbstractPotential)::AbstractFloat

Retrieve the radius cutoff for the provided potential.
This is the cutoff used for neighbor list calculations.
(i.e. Any pairs beyond this cutoff will be ignored.)
Defaults to `Inf` if a potential type does not implement a custom method.
"""
function get_rcutoff end
"""
    get_species(p::AbstractPotential)::Union{Tuple,Missing}

Retrieve the species to be included in an interaction
(pairs including a species not in the list are ignored).
A value of `missing` indicates that all species should be included,
which is the default behavior if a potential type does not implement a custom method.
"""
function get_species end
"""
    get_rpi(ace_params::ACE)::ACE1.RPI

Retrive the underlying RPI parameter type from ACE1. This is a convience function for exporting ACE parameters to file.
"""
function get_rpi end
"""
    length(params::BasisParameters) :: Int 

Retrieve the length of the descriptor vector, the number of descriptors used in the basis potential.
"""
function length end