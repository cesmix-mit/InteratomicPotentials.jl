################################################################################
# InteratomicPotentials API default generic implmentations
################################################################################

energy_and_force(A::AbstractSystem, p::ArbitraryPotential) = (; e = potential_energy(A, p), f = force(A, p))

virial(A::AbstractSystem, p::ArbitraryPotential) = sum(virial_stress(A, p)[1:3])
