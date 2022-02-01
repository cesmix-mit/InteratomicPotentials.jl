################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################

potential_energy(A::AbstractSystem, p::ArbitraryPotential) = energy_and_force(A, p).e
force(A::AbstractSystem, p::ArbitraryPotential) = energy_and_force(A, p).f
virial(A::AbstractSystem, p::ArbitraryPotential) = sum(virial_stress(A, p)[1:3])
