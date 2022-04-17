################################################################################
# InteratomicPotentials API implmentations for trainable potentials
################################################################################

function get_parameters(tp::TrainablePotential{P}) where {names,P<:Parameter{names}}
    P((getproperty(tp, name) for name ∈ names))
end
function set_parameters(tp::TP, p::P) where {P,TP<:TrainablePotential{P}}
    TP(; (name => getproperty(tp, name) for name ∈ propertynames(tp))..., p...)
end

function serialize_parameters(tp::TrainablePotential)
    collect(get_parameters(tp))
end
function deserialize_parameters(tp::TrainablePotential{P}, p::AbstractVector) where {P}
    set_parameters(tp, P(p))
end

function get_hyperparameters(tp::TrainablePotential{P,HP}) where {P,names,HP<:Parameter{names}}
    HP((getproperty(tp, name) for name ∈ names))
end
function set_hyperparameters(tp::TP, p::HP) where {P,HP,TP<:TrainablePotential{P,HP}}
    TP(; (name => getproperty(tp, name) for name ∈ propertynames(tp))..., p...)
end

function serialize_hyperparameters(tp::TrainablePotential)
    collect(get_hyperparameters(tp))
end
function deserialize_hyperparameters(tp::TrainablePotential{P,HP}, p::AbstractVector) where {P,HP}
    set_hyperparameters(tp, HP(p))
end