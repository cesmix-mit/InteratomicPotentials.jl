# The abstract types provided by InteratomicPotentials.jl

export AbstractPotential, NonTrainablePotential, TrainablePotential, EmpiricalPotential, MixedPotential

"""
    AbstractPotential

The abstract supertype of all interatomic potentials.
"""
abstract type AbstractPotential end
include("types/abstract_potential.jl")

"""
    NonTrainablePotential <: AbstractPotential

Abstract type for potentials that aren't trainable.
"""
abstract type NonTrainablePotential <: AbstractPotential end
include("types/non_trainable_potential.jl")

"""
    TrainablePotential{P<:NamedTuple,HP<:NamedTuple} <: AbstractPotential

Abstract type for potentials that are trainable.
`P` is a `NamedTuple` of parameter names and `HP`` is a `NamedTuple`` of hyperparameter names.
"""
abstract type TrainablePotential{P<:NamedTuple,HP<:NamedTuple} <: AbstractPotential end
include("types/trainable_potential.jl")

"""
    EmpiricalPotential{P<:NamedTuple,HP<:NamedTuple} <: TrainablePotential{P,HP}

# TODO (Dallas): describe what an empirical potential is
`P` is a `NamedTuple` of parameter names and `HP`` is a `NamedTuple`` of hyperparameter names.
"""
abstract type EmpiricalPotential{P<:NamedTuple,HP<:NamedTuple} <: TrainablePotential{P,HP} end
include("types/empirical_potential.jl")

"""
    MixedPotential <: AbstractPotential

Abstract type for potentials that are the combination of multiple sub-potentials.
"""
abstract type MixedPotential <: AbstractPotential end
include("types/mixed_potential.jl")
