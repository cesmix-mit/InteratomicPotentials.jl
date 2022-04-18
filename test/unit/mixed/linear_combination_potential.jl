@testset "Mixed Potential Unit Tests" begin
    atoms = [
        :Ar => (@SVector [0.0, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.5, 0.0, 0.0])u"bohr",
        :H => (@SVector [0.0, 0.5, 0.0])u"bohr",
        :H => (@SVector [0.5, 0.5, 0.0])u"bohr"
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    p1 = MockEmpiricalPotential(1.0, 1.0, (:Ar,))
    p2 = MockAbstractPotential(2.0)
    mixed_potentials = [
        +p1,
        -p1,
        p1 + p2,
        p1 - p2,
        2.0 * p1,
        p2 / 2.0,
        2.0 * p1 - 1.0 * p2,
        3.0 * (p1 + p2)
    ]

    true_rcutoffs = [
        1.0,
        1.0,
        Inf,
        Inf,
        1.0,
        Inf,
        Inf,
        Inf
    ]
    true_energies = [
        potential_energy(system, p1),
        -potential_energy(system, p1),
        potential_energy(system, p1) + potential_energy(system, p2),
        potential_energy(system, p1) - potential_energy(system, p2),
        2 * potential_energy(system, p1),
        potential_energy(system, p2) / 2,
        2 * potential_energy(system, p1) - potential_energy(system, p2),
        3 * (potential_energy(system, p1) + potential_energy(system, p2))
    ]
    true_forces = [
        force(system, p1),
        -force(system, p1),
        force(system, p1) + force(system, p2),
        force(system, p1) - force(system, p2),
        2 * force(system, p1),
        force(system, p2) / 2,
        2 * force(system, p1) - force(system, p2),
        3 * (force(system, p1) + force(system, p2))
    ]
    true_virials = [
        virial(system, p1),
        -virial(system, p1),
        virial(system, p1) + virial(system, p2),
        virial(system, p1) - virial(system, p2),
        2 * virial(system, p1),
        virial(system, p2) / 2,
        2 * virial(system, p1) - virial(system, p2),
        3 * (virial(system, p1) + virial(system, p2))
    ]
    true_virial_stresses = [
        virial_stress(system, p1),
        -virial_stress(system, p1),
        virial_stress(system, p1) + virial_stress(system, p2),
        virial_stress(system, p1) - virial_stress(system, p2),
        2 * virial_stress(system, p1),
        virial_stress(system, p2) / 2,
        2 * virial_stress(system, p1) - virial_stress(system, p2),
        3 * (virial_stress(system, p1) + virial_stress(system, p2))
    ]

    for (p, trc, te, tf, tv, tvs) in zip(mixed_potentials, true_rcutoffs, true_energies, true_forces, true_virials, true_virial_stresses)
        @test p isa MixedPotential

        @test get_rcutoff(p) == trc
        @test ismissing(get_species(p))

        @test energy_and_force(system, p).e == te
        @test energy_and_force(system, p).f == tf

        @test potential_energy(system, p) == te
        @test force(system, p) == tf
        @test virial(system, p) == tv
        @test virial_stress(system, p) == tvs
    end
end
