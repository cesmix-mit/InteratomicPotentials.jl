using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials
using Flux

# Define atomic system #########################################################

r0 = 2.0
position0 = @SVector [1.0, 0.0, 0.0]
element = :Ar
atom1 = Atom(element, position0 * u"Å")
atom2 = Atom(element, 0.25 * position0 * u"Å")
atom3 = Atom(element, 0.5 * position0 * u"Å")
atom4 = Atom(element, 0.75 * position0 * u"Å")
atom5 = Atom(element, 1.1 * position0 * u"Å")
atom6 = Atom(element, 2.0 * position0 * u"Å")
atoms = [atom1, atom2, atom3, atom4, atom5, atom6]
box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
bcs = [Periodic(), Periodic(), Periodic()]
system = FlexibleSystem(atoms, box * u"Å", bcs)

# Compose Neural ACE ###########################################################

# Define basis system
basis = ACE( species = [:Ar], body_order = 3, polynomial_degree = 3, 
             wL = 1.0, csp = 1.0, r0 = 0.4, rcutoff = 2.0)
n_desc = 5

# Define neural network besis potential
nns = Dict()
species = [element]
for s in species
    nns[s] = Chain( Dense(n_desc,128,σ; init = Flux.glorot_uniform(gain=-10)),
                    Dense(128,128,σ; init = Flux.glorot_uniform(gain=-10)),
                    Dense(128,1; init = Flux.glorot_uniform(gain=-10), bias = false))
end
nnbp = NNBasisPotential(nns, basis)
@test isa(nnbp, NNBasisPotential)

# Compute potential energy
@test isa(potential_energy(system, nnbp), Number)

# Compose Neural SNAP ##########################################################

# Define basis system
num_elements = 1
num_atoms = 3
twojmax = 4
rcutfac = 1.2
radii = [1.5]
rcut0 = 0.989
weight = [1.0]
chem_flag = false
bzero_flag = false
bnorm_flag = false
basis = SNAP(num_atoms, twojmax, [:Ar], rcutfac, 0.00, rcut0, radii, weight,
             chem_flag, bzero_flag, bnorm_flag)

n_desc = length(snap)

# Define neural network besis potential
nns = Dict()
species = [element]
for s in species
    nns[s] = Chain( Dense(n_desc,128,σ; init = Flux.glorot_uniform(gain=-10)),
                    Dense(128,128,σ; init = Flux.glorot_uniform(gain=-10)),
                    Dense(128,1; init = Flux.glorot_uniform(gain=-10), bias = false))
end
nnbp = NNBasisPotential(nns, basis)
@test isa(nnbp, NNBasisPotential)

# Compute potential energy
@test isa(potential_energy(system, nnbp), Number)

