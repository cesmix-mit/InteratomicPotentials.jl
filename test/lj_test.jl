using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position = @SVector [1.0, 0.0, 0.0] 
element  = :Ar
atom1     = StaticAtom(position * 1u"Å", element)

ϵ = 1.0
σ = 0.25
rcutoff  = 2.0
lj       = LennardJones(ϵ, σ, rcutoff)

@test isa(atom1, StaticAtom)
@test isa(lj, EmpiricalPotential)

atom2    = StaticAtom(0.25*position *1u"Å", element)
atom3    = StaticAtom(0.5*position *1u"Å", element)
atom4    = StaticAtom(0.75*position *1u"Å", element)
atom5    = StaticAtom(1.5*position *1u"Å", element)
atom6    = StaticAtom(2.0*position *1u"Å", element)
box = [[1., 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

system   = FlexibleSystem(box * 1u"Å", [DirichletZero(), Periodic(), Periodic()], [atom1, atom2, atom3, atom4, atom5, atom6])
n = length(system.particles)
nnlist   = neighborlist(system, rcutoff)
@test isa(nnlist, NeighborList)
println("NeighborList")
show(stdout, "text/plain", nnlist)
println("\n")

@test isa(potential_energy(system, lj), AbstractFloat)
println("Energy")
show(stdout, "text/plain", potential_energy(system, lj))
println("\n")
@test isa(force(system, lj), SVector{n, <:SVector{3, <:AbstractFloat}})

println("Force")
show(stdout, "text/plain", force(system, lj))
println("\n")

@test isa(virial(system, lj), AbstractFloat)
@test isa(virial_stress(system, lj), SVector{6, <:AbstractFloat})

println("Virial Stress")
show(stdout, "text/plain", virial_stress(system, lj))
println("\n")
