@testset "Mixed Potential Unit Tests" begin
    atoms = [
        :Ar => (@SVector [0.0, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.5, 0.0, 0.0])u"bohr",
        :H => (@SVector [0.0, 0.5, 0.0])u"bohr",
        :H => (@SVector [0.5, 0.5, 0.0])u"bohr"
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    ϵ = 1.0
    σ = 0.25
    rcutoff = 0.6
    lj1 = LennardJones(ϵ, σ, rcutoff, [:Ar, :H])
    lj2 = LennardJones(2ϵ, 2σ, rcutoff, [:Ar, :H])
    mixed_potentials = [
        +lj1,
        -lj1,
        lj1 + lj2,
        lj1 - lj2,
        2.0 * lj1,
        lj2 / 2.0,
        2.0 * lj1 - 1.0 * lj2,
        3.0 * (lj1 + lj2)
    ]

    true_energies = [
        potential_energy(system, lj1),
        -potential_energy(system, lj1),
        potential_energy(system, lj1) + potential_energy(system, lj2),
        potential_energy(system, lj1) - potential_energy(system, lj2),
        2.0 * potential_energy(system, lj1),
        potential_energy(system, lj2) / 2.0,
        2.0 * potential_energy(system, lj1) - 1.0 * potential_energy(system, lj2),
        3.0 * potential_energy(system, lj1) + 3.0 * potential_energy(system, lj2)
    ]
    true_forces = [
        force(system, lj1),
        -force(system, lj1),
        force(system, lj1) + force(system, lj2),
        force(system, lj1) - force(system, lj2),
        2.0 * force(system, lj1),
        force(system, lj2) / 2.0,
        2.0 * force(system, lj1) - 1.0 * force(system, lj2),
        3.0 * force(system, lj1) + 3.0 * force(system, lj2)
    ]
    true_virials = [
        virial(system, lj1),
        -virial(system, lj1),
        virial(system, lj1) + virial(system, lj2),
        virial(system, lj1) - virial(system, lj2),
        2.0 * virial(system, lj1),
        virial(system, lj2) / 2.0,
        2.0 * virial(system, lj1) - 1.0 * virial(system, lj2),
        3.0 * virial(system, lj1) + 3.0 * virial(system, lj2)
    ]
    true_virial_stresses = [
        virial_stress(system, lj1),
        -virial_stress(system, lj1),
        virial_stress(system, lj1) + virial_stress(system, lj2),
        virial_stress(system, lj1) - virial_stress(system, lj2),
        2.0 * virial_stress(system, lj1),
        virial_stress(system, lj2) / 2.0,
        2.0 * virial_stress(system, lj1) - 1.0 * virial_stress(system, lj2),
        3.0 * virial_stress(system, lj1) + 3.0 * virial_stress(system, lj2)
    ]

    for (result_lj, te, tf, tv, tvs) in zip(mixed_potentials, true_energies, true_forces, true_virials, true_virial_stresses)
        @test result_lj isa MixedPotential
        @test potential_energy(system, result_lj) isa ENERGY_TYPE
        @test potential_energy(system, result_lj) == te
        @test force(system, result_lj) isa AbstractVector{SVector{3,FORCE_TYPE}}
        @test force(system, result_lj) == tf
        @test virial(system, result_lj) isa ENERGY_TYPE
        @test virial(system, result_lj) == tv
        @test virial_stress(system, result_lj) isa SVector{6,ENERGY_TYPE}
        @test virial_stress(system, result_lj) == tvs

        e_f = energy_and_force(system, result_lj)
        @test e_f.e == te
        @test e_f.f == tf
    end
end
