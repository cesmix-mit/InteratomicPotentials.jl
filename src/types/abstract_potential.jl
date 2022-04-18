################################################################################
# InteratomicPotentials API default implementations
################################################################################

potential_energy(A::AbstractSystem, p::AbstractPotential) = energy_and_force(A, p).e
force(A::AbstractSystem, p::AbstractPotential) = energy_and_force(A, p).f
virial(A::AbstractSystem, p::AbstractPotential) = sum(virial_stress(A, p)[1:3])
