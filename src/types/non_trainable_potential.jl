################################################################################
# InteratomicPotentials API default implmentations for non-trainable potentials
################################################################################

get_parameters(::NonTrainablePotential) = Parameter{()}()
set_parameters(p::NonTrainablePotential, ::Parameter{()}) = p

serialize_parameters(::NonTrainablePotential) = []
deserialize_parameters(::NonTrainablePotential, ::AbstractVector) = p

get_hyperparameters(::NonTrainablePotential) = Parameter{()}()
set_hyperparameters(::NonTrainablePotential, ::Parameter{()}) = p

serialize_hyperparameters(::NonTrainablePotential) = []
deserialize_hyperparameters(::NonTrainablePotential, ::AbstractVector) = p
