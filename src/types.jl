# The abstract types provided by InteratomicPotentials.jl
import InteratomicPotentials: AbstractPotential
export BasisPotential, BasisParameters
"""
    BasisPotential <: AbstractPotential

Abstract type for potentials that are defined by a basis of polynomials.
"""
abstract type BasisPotential <: AbstractPotential end
"""
    BasisParameters 

Abstract type for the parameters that define the generation of descriptors for a potential.
"""
abstract type BasisParameters end
include("types/basis_potentials.jl")
