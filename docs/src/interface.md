# InteratomicBasisPotentials Interface

# Instantiating a Built-In Interatomic Potentials


Below is an example for building an Atomic Cluster Expansion (ACE) potential:
```julia
body_order = 2 
polynomial_degree = 8
r_inner_cutoff = 0.5 * u"Å"
r_cutoff = 4.0 * u"Å"
species = [:Ar, ]
wL      = 1.0 
csp     = 1.0
ace_params = RPIParams(species, body_order, polynomial_degree, wL, csp, r_inner_cutoff, r_cutoff)
```

