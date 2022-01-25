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
atom5 = Atom(element, 1.1 * position * u"Å")
atom6 = Atom(element, 2.0 * position * u"Å")
# atoms = [atom1, atom2, atom3, atom4, atom5, atom6]
atoms = [atom3, atom4, atom5]
box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
bcs = [Periodic(), Periodic(), Periodic()]
system = FlexibleSystem(atoms, box * u"Å", bcs)
n = length(system)
nnlist = neighborlist(system, rcutoff)
println(nnlist)
@test isa(nnlist, NeighborList)

true_r = [ [ SVector{3}([0.25, 0.0, 0.0]), SVector{3}([-.40, 0.0, 0.0]) ], 
            [ SVector{3}([-0.25, 0.0, 0.0]), SVector{3}([0.35, 0.0, 0.0]) ],
            [ SVector{3}([0.4, 0.0, 0.0]), SVector{3}([-0.35, 0.0, 0.0]) ] ]
nnlist_r = nnlist.r

for (ri, ri_true) in zip(nnlist_r, true_r)
    for (rij, rij_true) in zip(ri, ri_true)
        @test isapprox(sum(rij + rij_true), 0.0, atol = eps(2.0))
    end
end

@test isa(potential_energy(system, lj), AbstractFloat)
@test isa(force(system, lj), AbstractVector{<:SVector{3,<:AbstractFloat}})

@test isa(virial(system, lj), AbstractFloat)
@test isa(virial_stress(system, lj), SVector{6,<:AbstractFloat})


lj2 = LennardJones(2*ϵ, 2*σ, rcutoff)

sum_lj = lj+lj2
mixed_potentials = [lj+lj2, lj-lj2, 2.0*lj, -2.0*lj, 2.0*lj - 1.0*lj2, lj2 / 2.0]
true_energies = [potential_energy(system, lj) + potential_energy(system, lj2), 
                    potential_energy(system, lj) - potential_energy(system, lj2),
                    2.0*potential_energy(system, lj), 
                    -2.0*potential_energy(system, lj),
                    2.0*potential_energy(system, lj) - 1.0potential_energy(system, lj2), 
                    potential_energy(system, lj2)/2.0]

true_forces = [force(system, lj) + force(system, lj2), 
                    force(system, lj) - force(system, lj2),
                    2.0*force(system, lj), 
                    -2.0*force(system, lj),
                    2.0*force(system, lj) - 1.0force(system, lj2), 
                    force(system, lj2)/2.0]

true_virial = [virial(system, lj) + virial(system, lj2), 
                    virial(system, lj) - virial(system, lj2),
                    2.0*virial(system, lj), 
                    -2.0*virial(system, lj),
                    2.0*virial(system, lj) - 1.0virial(system, lj2), 
                    virial(system, lj2)/2.0]

for (result_lj, te, tf, tv) in zip(mixed_potentials, true_energies, true_forces, true_virial) 
    @test isa(result_lj, MixedPotential)
    @test isa(potential_energy(system, result_lj), AbstractFloat)
    @test isapprox(potential_energy(system, result_lj), te, rtol = 1e-6)
    @test isa(force(system, result_lj), AbstractVector{<:SVector{3,<:AbstractFloat}})
    @test isapprox(force(system, result_lj), tf, rtol = 1e-6)
    @test isa(virial(system, result_lj), AbstractFloat)
    @test isapprox(virial(system, result_lj), tv, rtol = 1e-6)
end
