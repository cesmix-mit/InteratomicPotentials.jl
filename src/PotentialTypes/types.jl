################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

abstract type ArbitraryPotential end
abstract type EmpiricalPotential <: ArbitraryPotential end
abstract type MixedPotential <: ArbitraryPotential end

################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
include("ArbitraryPotential/arbitrary_potential.jl")
export ArbitraryPotential

################################################################################
# Empirical Potentials
################################################################################
include("EmpiricalPotentials/empirical_potentials.jl")
export EmpiricalPotential

################################################################################
# Mixed Potentials
################################################################################
include("MixedPotentials/mixed_potentials.jl")
export MixedPotential

