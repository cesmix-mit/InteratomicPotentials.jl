# The abstract types provided by InteratomicPotentials.jl

export AbstractPotential, NonTrainablePotential, TrainablePotential, EmpiricalPotential, MixedPotential

################################################################################
# Abstract Supertype of all Potentials
################################################################################
abstract type AbstractPotential end
include("types/abstract_potential.jl")

################################################################################
# Non-Trainable Potentials
################################################################################
abstract type NonTrainablePotential <: AbstractPotential end
include("types/non_trainable_potential.jl")

################################################################################
# Trainable Potentials
################################################################################
abstract type TrainablePotential{P<:NamedTuple,HP<:NamedTuple} <: AbstractPotential end
include("types/trainable_potential.jl")

################################################################################
# Empirical Potentials
################################################################################
abstract type EmpiricalPotential{P<:NamedTuple,HP<:NamedTuple} <: TrainablePotential{P,HP} end
include("types/empirical_potential.jl")

################################################################################
# Mixed Potentials
################################################################################
abstract type MixedPotential <: AbstractPotential end
include("types/mixed_potential.jl")
