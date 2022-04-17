################################################################################
# InteratomicPotentials API default implementations for non-trainable potentials
################################################################################

get_parameters(::NonTrainablePotential) = NamedTuple{()}()
set_parameters(p::NonTrainablePotential, ::NamedTuple{()}) = p

serialize_parameters(::NonTrainablePotential) = []
deserialize_parameters(::NonTrainablePotential, ::AbstractVector) = p

get_hyperparameters(::NonTrainablePotential) = NamedTuple{()}()
set_hyperparameters(::NonTrainablePotential, ::NamedTuple{()}) = p

serialize_hyperparameters(::NonTrainablePotential) = []
deserialize_hyperparameters(::NonTrainablePotential, ::AbstractVector) = p
