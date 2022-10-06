# The InteratomicPotentials.jl API specification.

export energy_and_force, potential_energy, force, virial, virial_stress
export get_rcutoff, get_species
export get_parameters, set_parameters, serialize_parameters, deserialize_parameters
export get_hyperparameters, set_hyperparameters, serialize_hyperparameters, deserialize_hyperparameters
export compute_local_descriptors, compute_force_descriptors, compute_virial_descriptors, compute_all_descriptors
"""
    energy_and_force(s::AbstractSystem, p::AbstractPotential)::NamedTuple{(:e, :f), Tuple{Unitful.Energy,Vector{SVector{3, Unitful.Force}}}}

Calculate the unit-annotated potential energy of a system and the force acting on each particle in a system using the provided interatomic potential.
This combined function is offered because it is usually more efficient to calculate both properties simultaneously.
"""
function energy_and_force end
"""
    potential_energy(s::AbstractSystem, p::AbstractPotential)::Unitful.Energy

Calculate the unit-annotated potential energy of a system using the provided interatomic potential.
The default implementation uses the `:e` property of `energy_and_force(s,p)`.
"""
function potential_energy end
"""
    force(s::AbstractSystem, p::AbstractPotential)::Vector{SVector{3, Unitful.Force}}

Calculate the unit-annotated force acting on each particle in a system using the provided interatomic potential.
The default implementation uses the `:f` property of `energy_and_force(s,p)`.
"""
function force end
"""
    virial_stress(s::AbstractSystem, p::AbstractPotential)::SVector{6,Unitful.Energy}

Calculate the unit-annotated virial stress tensor of a system, officially calculated as the sum of radial-force outerproducts: `` \\sum r_{ij} \\bigotimes F_{ij}``, only returns the unique lower-diagonal components.
"""
function virial_stress end
"""
    virial(s::AbstractSystem, p::AbstractPotential)::Unitful.Energy

Calculate the unit-annotated virial of a system, officially calculated as the trace contraction of the sum of radial-force outerproducts: ``tr\\left( \\sum r_{ij} \\bigotimes F_{ij} \\right)``
"""
function virial end
"""
    get_rcutoff(p::AbstractPotential)::AbstractFloat

Retrieve the radius cutoff for the provided potential.
This is the cutoff used for neighbor list calculations.
(i.e. Any pairs beyond this cutoff will be ignored.)
Defaults to `Inf` if a potential type does not implement a custom method.
"""
function get_rcutoff end
"""
    get_species(p::AbstractPotential)::Union{Tuple,Missing}

Retrieve the species to be included in an interaction
(pairs including a species not in the list are ignored).
A value of `missing` indicates that all species should be included,
which is the default behavior if a potential type does not implement a custom method.
"""
function get_species end

"""
    get_parameters(p::AbstractPotential)::NamedTuple

Retrieve a `NamedTuple` with the trainable parameters of a potential.

For a `TrainablePotential{P,HP}`, the function automatically extracts the fields with names in the type parameter `P`.
For a `NonTrainablePotential`, the function returns an empty `NamedTuple`.
"""
function get_parameters end
"""
    set_parameters(p::P<:AbstractPotential, parameters::NamedTuple)::P

Generate a new potential of the same type as an exising potential but with modified trainable parameters.

For a `TrainablePotential{P,HP}`, the function copies over all fields and replaces the fields from the type parameter `P`.
For a `NonTrainablePotential`, the function returns the original potential, unchanged.
"""
function set_parameters end
"""
    serialize_parameters(p::AbstractPotential)::Vector{<:AbstractFloat}

Generate a `Vector` form of the trainable parameters of a potential.

For a `TrainablePotential{P,HP}`, the function automatically extracts the fields with names in the type parameter `P` and collects them into a `Vector`.
Note that this behavior should be overloaded for a potential type that has any non-scalar trainable parameters.
For a `NonTrainablePotential`, the function returns an empty `Vector`.
"""
function serialize_parameters end
"""
    deserialize_parameters(p::P<:AbstractPotential, parameters::AbstractVector{<:AbstractFloat})::P

Generate a new potential of the same type as an exising potential but with modified trainable parameters.
This function takes in the trainable parameters as a `Vector` in the same form as produced by `serialize_parameters(p)`.

For a `TrainablePotential{P,HP}`, the function converts the `Vector` of parameters into the type parameter `P` then calls the `set_parameters` implementation.
Note that this behavior should be overloaded for a potential type that has any non-scalar trainable parameters.
For a `NonTrainablePotential`, the function returns the original potential, unchanged.
"""
function deserialize_parameters end

"""
    get_hyperparameters(p::AbstractPotential)::NamedTuple

Retrieve a `NamedTuple` with the trainable hyperparameters of a potential.

For a `TrainablePotential{P,HP}`, the function automatically extracts the fields with names in the type parameter `HP`.
For a `NonTrainablePotential`, the function returns an empty `NamedTuple`.
"""
function get_hyperparameters end
"""
    set_hyperparameters(p::P<:AbstractPotential, hyperparameters::NamedTuple)::P

Generate a new potential of the same type as an exising potential but with modified trainable hyperparameters.

For a `TrainablePotential{P,HP}`, the function copies over all fields and replaces the fields from the type parameter `HP`.
For a `NonTrainablePotential`, the function returns the original potential, unchanged.
"""
function set_hyperparameters end
"""
    serialize_hyperparameters(p::AbstractPotential)::Vector{<:AbstractFloat}

Generate a `Vector` form of the trainable hyperparameters of a potential.

For a `TrainablePotential{P,HP}`, the function automatically extracts the fields with names in the type parameter `HP` and collects them into a `Vector`.
Note that this behavior should be overloaded for a potential type that has any non-scalar trainable hyperparameters.
For a `NonTrainablePotential`, the function returns an empty `Vector`.
"""
function serialize_hyperparameters end
"""
    deserialize_hyperparameters(p::P<:AbstractPotential, hyperparameters::AbstractVector{<:AbstractFloat})::P

Generate a new potential of the same type as an exising potential but with modified trainable hyperparameters.
This function takes in the trainable hyperparameters as a `Vector` in the same form as produced by `serialize_hyperparameters(p)`.

For a `TrainablePotential{P,HP}`, the function converts the `Vector` of hyperparameters into the type parameter `HP` then calls the `set_hyperparameters` implementation.
Note that this behavior should be overloaded for a potential type that has any non-scalar trainable hyperparameters.
For a `NonTrainablePotential`, the function returns the original potential, unchanged.
"""
function deserialize_hyperparameters end
"""
    get_local_descriptors(s::AbstractSystem, b::BasisSystem)::Vector{Vector}

Compute the local descriptors for an abstract system s using the basis system defined in b. These types are defined in this package, but the primary use is in the InteratomicBasisPotentials.jl package, where the individual basis systems are defined. 

For a system s with n_atom number of atoms this method get_local_descriptors returns a Vector{Vector{T<:Real}}, containing n_atom vectors each of dimension d as defined by the basis system.
"""
function compute_local_descriptors end 
"""
    compute_force_descriptors(s::AbstractSystem, b::BasisSystem)::Vector{Vector{Vector}}

Compute the force descriptors for an abstract system s using the basis system defined in b. These types are defined in this package, but the primary use is in the InteratomicBasisPotentials.jl package, where the individual basis systems are defined. 

For a system s with n_atom number of atoms this method compute_force_descriptors returns a Vector{Matrix{T<:Real}}}, containing n_atom vectors each containing a matrix with rows corresponding to x-, y-, and z-components of dimension d as defined by the basis system.
"""
function compute_force_descriptors end
"""
    compute_virial_descriptors(s::AbstractSystem, b::BasisSystem)::Vector{Vector{Vector}}

Compute the virial descriptors for an abstract system s using the basis system defined in b. These types are defined in this package, but the primary use is in the InteratomicBasisPotentials.jl package, where the individual basis systems are defined. 

For a system s with n_atom number of atoms this method compute_virial_descriptors returns a Matrix{T<:Real}, containing 6 rows corresponding to entries of the stress tensor and columns of dimension d as defined by the basis system.
"""
function compute_virial_descriptors end
"""
    compute_all_descriptors(s::AbstractSystem, b::BasisSystem)::Vector{Vector{Vector}}

Compute the local, force, and virial descriptors for an abstract system s using the basis system defined in b. These types are defined in this package, but the primary use is in the InteratomicBasisPotentials.jl package, where the individual basis systems are defined. 
"""
function compute_all_descriptors end

