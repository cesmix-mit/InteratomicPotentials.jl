@testset "Mixed Potential Unit Tests" begin
    atoms = [
        :Ar => (@SVector [0.0, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.5, 0.0, 0.0])u"bohr",
        :H => (@SVector [0.0, 0.5, 0.0])u"bohr",
        :H => (@SVector [0.5, 0.5, 0.0])u"bohr"
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    p1 = MockArbitraryPotential(1.0, (:Ar, :H))
    p2 = MockArbitraryPotential(2.0, (:Ar, :H))
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

    true_energies = [
        1.0,
        -1.0,
        3.0,
        -1.0,
        2.0,
        1.0,
        0.0,
        9.0
    ]u"hartree"
    true_forces = [
        SVector{3}.([[2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [2.0, 3.0, 4.0]]),
        SVector{3}.([[-2.0, -3.0, -4.0], [-2.0, -3.0, -4.0], [-2.0, -3.0, -4.0], [-2.0, -3.0, -4.0]]),
        SVector{3}.([[5.0, 7.0, 9.0], [5.0, 7.0, 9.0], [5.0, 7.0, 9.0], [5.0, 7.0, 9.0]]),
        SVector{3}.([[-1.0, -1.0, -1.0], [-1.0, -1.0, -1.0], [-1.0, -1.0, -1.0], [-1.0, -1.0, -1.0]]),
        SVector{3}.([[4.0, 6.0, 8.0], [4.0, 6.0, 8.0], [4.0, 6.0, 8.0], [4.0, 6.0, 8.0]]),
        SVector{3}.([[1.5, 2.0, 2.5], [1.5, 2.0, 2.5], [1.5, 2.0, 2.5], [1.5, 2.0, 2.5]]),
        SVector{3}.([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]),
        SVector{3}.([[15.0, 21.0, 27.0], [15.0, 21.0, 27.0], [15.0, 21.0, 27.0], [15.0, 21.0, 27.0]])
    ]u"hartree/bohr"
    true_virials = [
        -3.0,
        3.0,
        -3.0,
        -3.0,
        -6.0,
        0.0,
        -6.0,
        -9.0
    ]u"hartree"
    true_virial_stresses = SVector{6}.([
        [0.0, -1.0, -2.0, -3.0, -4.0, -5.0],
        [-0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [1.0, -1.0, -3.0, -5.0, -7.0, -9.0],
        [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        [0.0, -2.0, -4.0, -6.0, -8.0, -10.0],
        [0.5, 0.0, -0.5, -1.0, -1.5, -2.0],
        [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
        [3.0, -3.0, -9.0, -15.0, -21.0, -27.0]
    ])u"hartree"

    for (p, te, tf, tv, tvs) in zip(mixed_potentials, true_energies, true_forces, true_virials, true_virial_stresses)
        @test p isa MixedPotential

        @test potential_energy(system, p) == te
        @test force(system, p) == tf
        @test virial(system, p) == tv
        @test virial_stress(system, p) == tvs

        @test energy_and_force(system, p).e == te
        @test energy_and_force(system, p).f == tf
    end
end
