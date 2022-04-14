################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
abstract type ArbitraryPotential end
include("ArbitraryPotential/arbitrary_potential.jl")
export ArbitraryPotential

################################################################################
# Empirical Potentials
################################################################################
abstract type EmpiricalPotential <: ArbitraryPotential end
include("EmpiricalPotentials/empirical_potentials.jl")
export EmpiricalPotential

################################################################################
# Mixed Potentials
################################################################################
abstract type MixedPotential <: ArbitraryPotential end
include("MixedPotentials/mixed_potentials.jl")
export MixedPotential
