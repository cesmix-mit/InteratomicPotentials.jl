################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

abstract type ArbitraryPotential end
abstract type EmpiricalPotential <: ArbitraryPotential end
abstract type BasisPotential <: ArbitraryPotential end
abstract type MixedPotential <: ArbitraryPotential end

################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
include("ArbitraryPotential/arbitrary_potential.jl")

################################################################################
# Empirical Potentials
################################################################################
include("EmpiricalPotentials/empirical_potentials.jl")

################################################################################
# BasisPotentials
################################################################################
include("BasisPotentials/basis_potentials.jl")

################################################################################
# GaN
################################################################################
# include("GaN/gan.jl")
