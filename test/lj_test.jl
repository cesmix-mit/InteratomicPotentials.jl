using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position = @SVector [1.0, 0.0, 0.0]
element = :Ar
atom1 = Atom(element, position * u"Å")

ϵ = 1.0
σ = 0.25
rcutoff = 2.0
lj = LennardJones(ϵ, σ, rcutoff)

@test isa(atom1, Atom)
@test isa(lj, EmpiricalPotential)

atom2 = Atom(element, 0.25 * position * u"Å")
atom3 = Atom(element, 0.5 * position * u"Å")
atom4 = Atom(element, 0.75 * position * u"Å")
atom5 = Atom(element, 1.5 * position * u"Å")
atom6 = Atom(element, 2.0 * position * u"Å")
atoms = [atom1, atom2, atom3, atom4, atom5, atom6]
box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
bcs = [DirichletZero(), Periodic(), Periodic()]
system = FlexibleSystem(atoms, box * u"Å", bcs)
n = length(system)
nnlist = neighborlist(system, rcutoff)
@test isa(nnlist, NeighborList)

@test isa(potential_energy(system, lj), AbstractFloat)
@test isa(force(system, lj), AbstractVector{<:SVector{3,<:AbstractFloat}})

@test isa(virial(system, lj), AbstractFloat)
@test isa(virial_stress(system, lj), SVector{6,<:AbstractFloat})

@test potential_energy(system, lj) == InteratomicPotentials.potential_energy2(system, lj)
@test force(system, lj) == InteratomicPotentials.force2(system, lj)
@test virial(system, lj) == InteratomicPotentials.virial2(system, lj)
@test virial_stress(system, lj) == InteratomicPotentials.virial_stress2(system, lj)
