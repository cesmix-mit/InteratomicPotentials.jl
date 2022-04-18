################################################################################
# InteratomicPotentials API default implementations
################################################################################

get_rcutoff(::AbstractPotential) = Inf
get_species(::AbstractPotential) = missing

potential_energy(s::AbstractSystem, p::AbstractPotential) = energy_and_force(s, p).e
force(s::AbstractSystem, p::AbstractPotential) = energy_and_force(s, p).f
virial(s::AbstractSystem, p::AbstractPotential) = sum(virial_stress(s, p)[1:3])
