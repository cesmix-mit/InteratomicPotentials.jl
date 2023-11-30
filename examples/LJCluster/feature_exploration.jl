push!(Base.LOAD_PATH, "../../")

using InteratomicPotentials 
using AtomsBase
using UnitfulAtomic
using Unitful 
using StaticArrays
using LinearAlgebra 
using Statistics 

# Load data 
include("load_data.jl")
systems, energies, forces = load_data(;num_entries = 100);

# Split into training, testing
train_systems, train_energies, train_forces = systems[1:50], energies[1:50], forces[1:50];
test_systems, test_energies, test_forces = systems[51:end], energies[51:end], forces[51:end];

## Create ACE Basis (2body, 8 polynomial degree)
ace = ACE( species = [:Ar],
           body_order = 2,
           polynomial_degree = 8, 
           wL = 1.0,
           csp = 1.0,
           r0 = 1.0,
           rcutoff = 2.0)

## Calculate descriptors 
train_descriptors = [sum(compute_local_descriptors(sys, ace)) for sys in train_systems]
test_descriptors = [sum(compute_local_descriptors(sys, ace)) for sys in test_systems]
mean_d = mean([train_descriptors; test_descriptors])
cov_d = cov([train_descriptors; test_descriptors])
inv_cov_d = inv(cov_d)
mah_dist_train = [ sqrt( 0.5*(td - mean_d)' * inv_cov_d * (td - mean_d) ) for td in train_descriptors ]
mah_dist_test = [ sqrt( 0.5*(td - mean_d)' * inv_cov_d * (td - mean_d) ) for td in test_descriptors ]

# Create linear basis potential and estimate β
lb = LBasisPotential(ace)
lb.β .= hcat(train_descriptors...)' \ train_energies

# Calculate potential energy
potential_energy(test_systems[1], lb)

