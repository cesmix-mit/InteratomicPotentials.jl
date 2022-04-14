using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

ϵ = 1.0
σ = 0.25
rcutoff = 2.0
lj = LennardJones(ϵ * u"hartree", σ * u"bohr", rcutoff * u"bohr", [:Ar, :H])

@test lj isa EmpiricalPotential

position = (@SVector [1.0, 0.0, 0.0])u"bohr"
element = :Ar

atom1 = Atom(:Ar, position)
atom2 = Atom(:Ar, 0.25 * position)
atom3 = Atom(:Ar, 0.5 * position)
atom4 = Atom(:Ar, 0.75 * position)
atom5 = Atom(:H, 1.1 * position)
atom6 = Atom(:Ar, 2.0 * position)
atoms = [atom3, atom4, atom5]
box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
bcs = [Periodic(), Periodic(), Periodic()]
system = FlexibleSystem(atoms, box, bcs)

nnlist = neighborlist(system, rcutoff)

nnlist_r = nnlist.r
true_r = [
    [SVector{3}([0.25, 0.0, 0.0]), SVector{3}([-0.40, 0.0, 0.0])],
    [SVector{3}([0.35, 0.0, 0.0])],
    []
]

for (ri, ri_true) in zip(nnlist_r, true_r)
    for (rij, rij_true) in zip(ri, ri_true)
        @test isapprox(sum(rij + rij_true), 0.0, atol=eps(2.0))
    end
end

ENERGY_TYPE = typeof(0.0u"hartree")
FORCE_TYPE = typeof(0.0u"hartree/bohr")

@test potential_energy(system, lj) isa ENERGY_TYPE
@test force(system, lj) isa AbstractVector{SVector{3,FORCE_TYPE}}

@test virial(system, lj) isa ENERGY_TYPE
@test virial_stress(system, lj) isa SVector{6,ENERGY_TYPE}

lj2 = LennardJones(2 * ϵ, 2 * σ, rcutoff, [:Ar, :H])

sum_lj = lj + lj2
mixed_potentials = [lj + lj2, lj - lj2, 2.0 * lj, -2.0 * lj, 2.0 * lj - 1.0 * lj2, lj2 / 2.0]

true_energies = [
    potential_energy(system, lj) + potential_energy(system, lj2),
    potential_energy(system, lj) - potential_energy(system, lj2),
    2.0 * potential_energy(system, lj),
    -2.0 * potential_energy(system, lj),
    2.0 * potential_energy(system, lj) - 1.0 * potential_energy(system, lj2),
    potential_energy(system, lj2) / 2.0
]
true_forces = [
    force(system, lj) + force(system, lj2),
    force(system, lj) - force(system, lj2),
    2.0 * force(system, lj),
    -2.0 * force(system, lj),
    2.0 * force(system, lj) - 1.0 * force(system, lj2),
    force(system, lj2) / 2.0
]
true_virial = [
    virial(system, lj) + virial(system, lj2),
    virial(system, lj) - virial(system, lj2),
    2.0 * virial(system, lj),
    -2.0 * virial(system, lj),
    2.0 * virial(system, lj) - 1.0 * virial(system, lj2),
    virial(system, lj2) / 2.0
]

for (result_lj, te, tf, tv) in zip(mixed_potentials, true_energies, true_forces, true_virial)
    @test result_lj isa MixedPotential
    @test potential_energy(system, result_lj) isa ENERGY_TYPE
    @test potential_energy(system, result_lj) == te
    @test force(system, result_lj) isa AbstractVector{SVector{3,FORCE_TYPE}}
    @test force(system, result_lj) == tf
    @test virial(system, result_lj) isa ENERGY_TYPE
    @test virial(system, result_lj) == tv

    e_f = energy_and_force(system, result_lj)
    @test e_f.e == te
    @test e_f.f == tf
end
