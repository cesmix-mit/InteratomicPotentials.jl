using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position = @SVector [1.0, 1.0, 1.0] 
element  = :Ar
atom     = StaticAtom(position * 1u"Å", element)
lj       = LennardJones(1.0, 1.0)
@test isa(atom, StaticAtom)
@test isa(lj, EmpiricalPotential)
@test isa(potential_energy(atom, lj), AbstractFloat)
@test isa(force(atom, lj), SVector{3, <:AbstractFloat})
@test isa(virial(atom, lj), AbstractFloat)
@test isa(virial_stress(atom, lj), SVector{6, <:AbstractFloat})

atom2    = StaticAtom(0.0*position *1u"Å", element)
box = [[0.0, 1.0], [0.0, 1.0]]
system   = FlexibleSystem(box * 1u"Å", [Periodic(), Periodic()], [atom, atom2])

@test isa(potential_energy(system, lj), AbstractFloat)
@test isa(force(system, lj), SVector{2, <:SVector{3, <:AbstractFloat}})
@test isa(virial(system, lj), AbstractFloat)
@test isa(virial_stress(system, lj), SVector{6, <:AbstractFloat})