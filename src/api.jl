# The InteratomicPotentials.jl API specification.

export energy_and_force, potential_energy, force, virial, virial_stress
export get_rcutoff, get_species
export get_parameters, set_parameters, serialize_parameters, deserialize_parameters
export get_hyperparameters, set_hyperparameters, serialize_hyperparameters, deserialize_hyperparameters

function energy_and_force end
function potential_energy end
function force end
function virial end
function virial_stress end

function get_rcutoff end
function get_species end

function get_parameters end
function set_parameters end
function serialize_parameters end
function deserialize_parameters end

function get_hyperparameters end
function set_hyperparameters end
function serialize_hyperparameters end
function deserialize_hyperparameters end
