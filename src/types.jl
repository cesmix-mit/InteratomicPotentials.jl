# The abstract types provided by InteratomicPotentials.jl

export ArbitraryPotential, NonTrainablePotential, TrainablePotential, EmpiricalPotential, MixedPotential

################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################
abstract type ArbitraryPotential end
include("types/arbitrary_potential.jl")

################################################################################
# Non-Trainable Potentials
################################################################################
abstract type NonTrainablePotential <: ArbitraryPotential end
include("types/non_trainable_potential.jl")

################################################################################
# Trainable Potentials
################################################################################
abstract type TrainablePotential{P<:Parameter,HP<:Parameter} <: ArbitraryPotential end
include("types/trainable_potential.jl")

################################################################################
# Empirical Potentials
################################################################################
abstract type EmpiricalPotential{P<:Parameter,HP<:Parameter} <: TrainablePotential{P,HP} end
include("types/empirical_potential.jl")

################################################################################
# Mixed Potentials
################################################################################
abstract type MixedPotential <: ArbitraryPotential end
include("types/mixed_potential.jl")
