import Base.NamedTuple as Parameter
export Parameter

################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
abstract type ArbitraryPotential end
include("ArbitraryPotentials/arbitrary_potential.jl")
export ArbitraryPotential

################################################################################
# Empirical Potentials
################################################################################
abstract type EmpiricalPotential <: ArbitraryPotential end
include("EmpiricalPotentials/empirical_potential.jl")
export EmpiricalPotential

################################################################################
# Mixed Potentials
################################################################################
abstract type MixedPotential <: ArbitraryPotential end
include("MixedPotentials/mixed_potential.jl")
export MixedPotential
