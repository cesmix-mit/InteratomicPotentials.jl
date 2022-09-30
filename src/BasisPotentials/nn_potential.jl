"""
NNBasisPotential

Definition of the neural network basis potential composed type.
"""
struct NNBasisPotential <: BasisPotential{NamedTuple{(:params, )}, NamedTuple{()}}
    nn :: Union{Flux.Chain, Flux.Dense}
    basis :: BasisSystem
end

"""
potential_energy(A::AbstractSystem, p::NNBasisPotential)

`A`: atomic system or atomic configuration.
`p`: neural network basis potential.

Returns the potential energy of a system using a neural network basis potential.
"""
function potential_energy(A::AbstractSystem, p::NNBasisPotential)
    sum(p.nn(sum(get_local_descriptors(A, p.basis))))
end


"""
force(A::AbstractSystem, p::NNBasisPotential)

`A`: atomic system or atomic configuration.
`p`: neural network basis potential.

Returns the force of a system using a neural network basis potential.
"""
function force(A::AbstractSystem, p::NNBasisPotential)
    b = sum(get_local_descriptors(A, p.basis))
    dnndb = sum(gradient(p.nn, b))
    dbdr = get_force_descriptors(A, p.basis)
    return [- dbdr[atom] * dnndb  
         for atom in 1:length(dbdr)]
end


get_parameters(p::NNBasisPotential) = Flux.params(p.nn)

# """
# potential_energy(b::Vector, ps::Vector, re)

# `b`: energy descriptors of a system.
# `ps`: neural network parameters. See Flux.destructure.
# `re`: neural network restructure. See Flux.destructure.

# Returns the potential energy of a system using a neural network basis potential.
# """
# function potential_energy(b::Vector, ps::Vector, re)
# return sum(re(ps)(b))
# end


# TODO: create issue.
#Note: calculating the gradient of the loss function requires in turn
#calculating the gradient of the energy. That is, calculating the gradient of
#a function that calculates another gradient.
#So far I have not found a clean way to do this using the abstractions 
#provided by Flux, which in turn is integrated with Zygote. 
#Links related to this issue:
#https://discourse.julialang.org/t/how-to-add-norm-of-gradient-to-a-loss-function/69873/16
#https://discourse.julialang.org/t/issue-with-zygote-over-forwarddiff-derivative/70824
#https://github.com/FluxML/Zygote.jl/issues/953#issuecomment-841882071
#
#To solve this for the moment I am calculating one of the gradients analytically.
#To do this I had to use Flux.destructure, which I think makes the code slower
#because of the copies it creates.


# """
# grad_mlp(nn_params, x0)

# `nn_params`: neural network parameters.
# `x0`: first layer input (energy descriptors).

# Returns the analytical derivative of a feedforward neural network.
# """
# function grad_mlp(nn_params, x0)
# dsdy(x) = x>0 ? 1 : 0 # Flux.σ(x) * (1 - Flux.σ(x)) 
# prod = 1; x = x0
# n_layers = length(nn_params) ÷ 2
# for i in 1:2:2(n_layers-1)-1  # i : 1, 3
#     y = nn_params[i] * x + nn_params[i+1]
#     x =  Flux.relu.(y) # Flux.σ.(y)
#     prod = dsdy.(y) .* nn_params[i] * prod
# end
# i = 2(n_layers)-1 
# prod = nn_params[i] * prod
# return prod
# end


# """
# force(b::Vector, dbdr::Vector, p::NNBasisPotential)

# `b`: energy descriptors of a system.
# `dbdr`: force descriptors of a system.
# `p`: neural network basis potential.

# Returns the force of a system using a neural network basis potential.
# """
# function force(b::Vector, dbdr::Vector, p::NNBasisPotential)
# dnndb = grad_mlp(p.nn_params, b)
# return dnndb ⋅ dbdr
# end


# """
# force(b::Vector, dbdr::Vector, ps::Vector, re)

# `b`: energy descriptors of a system.
# `dbdr`: force descriptors of a system.
# `ps`: neural network parameters. See Flux.destructure.
# `re`: neural network restructure. See Flux.destructure.

# Returns the force of a system using a neural network basis potential.
# """
# function force(b::Vector, dbdr::Vector, ps::Vector, re)
# nn_params = Flux.params(re(ps))
# dnndb = grad_mlp(nn_params, b)
# return dnndb ⋅ dbdr
# end

# Computing the force using ForwardDiff
#function force(b::Vector, dbdr::Vector, p::NNBasisPotential)
#    dnndb = ForwardDiff.gradient(x -> sum(p.nn(x)), b)
#    return dnndb ⋅ dbdr
#end

#function force(b::Vector, dbdr::Vector, p::NNBasisPotential)
#    dnndb = gradient(x -> sum(p.nn(x)), b)[1]
#    return dnndb ⋅ dbdr
#end

# Computing the force using ForwardDiff and destructure
#function force(b::Vector, dbdr::Vector, ps::Vector, re)
#    dnndb = ForwardDiff.gradient(x -> sum(re(ps)(x)), b)
#    return dnndb ⋅ dbdr
#end

# Computing the force using pullback
#function force(b::Vector, dbdr::Vector, p::NNBasisPotential)
#    y, pullback = Zygote.pullback(p.nn, b)
#    dnndb = pullback(ones(size(y)))[1]
#    return dnndb ⋅ dbdr
#end

#function force(b::Vector, dbdr::Vector,  ps::Vector, re)
#    y, pullback = Zygote.pullback(re(ps), b)
#    dnndb = pullback(ones(size(y)))[1]
#    return dnndb ⋅ dbdr
#end

