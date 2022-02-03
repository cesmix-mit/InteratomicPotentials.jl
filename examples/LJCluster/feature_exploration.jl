using InteratomicPotentials 
using AtomsBase
using UnitfulAtomic
using Unitful 
using StaticArrays
using LinearAlgebra 
using Statistics 

# Load data 
include("LJCluster/load_data.jl")
systems, energies, forces = load_data(;num_entries = 100);

# Split into training, testing
train_systems, train_energies, train_forces = systems[1:50], energies[1:50], forces[1:50];
test_systems, test_energies, test_forces = systems[51:end], energies[51:end], forces[51:end];

## Create RPI Basis (2body, 8 polynomial degree)
n_body = 2
max_deg = 8
r0 = 1.0
rcutoff = 4.0
wL = 1.0
csp = 1.0
rpi_params = RPIParams([:Ar], n_body, max_deg, wL, csp, r0, rcutoff)

## Calculate descriptors 
train_descriptors = [evaluate_basis(sys, rpi_params) for sys in train_systems]
test_descriptors = [evaluate_basis(sys, rpi_params) for sys in test_systems]
mean_d = mean(descriptors)
cov_d = cov(descriptors)
inv_cov_d = inv(cov_d)
mah_dist_train = [ sqrt( 0.5*(td - mean_d)' * inv_cov_d * (td - mean_d) ) for td in train_descriptors ]
mah_dist_test = [ sqrt( 0.5*(td - mean_d)' * inv_cov_d * (td - mean_d) ) for td in test_descriptors ]

# Estimate β
β = hcat(train_descriptors...)' \ train_energies
rpi = RPI(β, rpi_params)

potential_energy(test_systems[1], rpi)