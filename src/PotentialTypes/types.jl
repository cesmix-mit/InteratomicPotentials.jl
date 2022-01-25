################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

abstract type ArbitraryPotential end
abstract type ElementalPotential <: ArbitraryPotential end # Single potentials
abstract type EmpiricalPotential <: ElementalPotential end
abstract type BasisPotential <: ElementalPotential end

abstract type MixedPotential <: ArbitraryPotential end # Combined potentials
############################### Empirical Potentials ################################################
include("EmpiricalPotentials/empirical_potentials.jl")

################################ BasisPotentials ##############################################################
include("BasisPotentials/basis_potentials.jl")

################################ MixedPotential ###############################################################
include("MixedPotentials/mixed_potentials.jl")
