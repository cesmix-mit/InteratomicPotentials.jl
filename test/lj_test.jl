using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position = @SVector [1.0, 1.0, 1.0] 
element  = ChemicalElement(:Ar)
atom     = SimpleAtom(position * 1u"Å", element)
lj       = LennardJones(1.0, 1.0)
@test isa(atom, SimpleAtom)
@test isa(lj, EmpiricalPotential)
@test isa(potential_energy(atom, lj), AbstractFloat)
@test isa(force(atom, lj), SVector{3, <:AbstractFloat})
@test isa(virial(atom, lj), AbstractFloat)
@test isa(virial_stress(atom, lj), SVector{6, <:AbstractFloat})

atom2    = SimpleAtom(0.0*position *1u"Å", element)
box = SVector(SVector([0.0, 1.0], [0.0, 1.0]))
system   = SimpleSystem(box * 1u"Å", SVector{2}([Periodic(), Periodic()]), [atom, atom2])

@test isa(potential_energy(system, lj), AbstractFloat)
print(typeof(force(system, lj)))
@test isa(force(system, lj), SVector{2, <:SVector{3, <:AbstractFloat}})
@test isa(virial(system, lj), AbstractFloat)
@test isa(virial_stress(system, lj), SVector{6, <:AbstractFloat})