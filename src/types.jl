# The abstract types provided by InteratomicPotentials.jl

export ArbitraryPotential, EmpiricalPotential, MixedPotential

################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
abstract type ArbitraryPotential end
include("types/arbitrary_potential.jl")

################################################################################
# Empirical Potentials
################################################################################
abstract type EmpiricalPotential <: ArbitraryPotential end
include("types/empirical_potential.jl")

################################################################################
# Mixed Potentials
################################################################################
abstract type MixedPotential <: ArbitraryPotential end
include("types/mixed_potential.jl")
