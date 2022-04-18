@testset "Abstract Potential Default Implmentation Unit Tests" begin
    atoms = [
        :Ar => (@SVector [0.0, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.5, 0.0, 0.0])u"bohr",
        :H => (@SVector [0.0, 0.5, 0.0])u"bohr",
        :H => (@SVector [0.5, 0.5, 0.0])u"bohr"
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    p = MockAbstractPotential(0.0)

    @test get_rcutoff(p) == Inf
    @test ismissing(get_species(p))

    @test potential_energy(system, p) == 0.0u"hartree"
    @test force(system, p) == fill((@SVector[1.0, 2.0, 3.0])u"hartree/bohr", 4)
    @test virial(system, p) == -6u"hartree"
end
