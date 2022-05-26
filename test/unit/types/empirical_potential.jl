@testset "Empirical Potential Default Implmentation Unit Tests" begin
    atoms = [
        :Ar => (@SVector [0.0, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.5, 0.0, 0.0])u"bohr",
        :Ar => (@SVector [0.0, 0.5, 0.0])u"bohr",
        :H => (@SVector [0.5, 0.5, 0.0])u"bohr"
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    p = MockEmpiricalPotential(1.0, 0.5, (:Ar,))

    @test potential_energy(@SVector[1.0, 1.0, 1.0], p) ≈ 3.0
    @test force(@SVector[1.0, 1.0, 1.0], p) ≈ @SVector[3.0, 3.0, 3.0]

    @test energy_and_force(system, p) isa NamedTuple{(:e, :f),<:Tuple{<:Unitful.Energy,<:AbstractVector{<:SVector{3,<:Unitful.Force}}}}
    @test energy_and_force(system, p).e == 0.5u"hartree"
    @test energy_and_force(system, p).f == [
        (@SVector[-0.125, -0.125, 0.0])u"hartree/bohr",
        (@SVector[0.125, 0.0, 0.0])u"hartree/bohr",
        (@SVector[0.0, 0.125, 0.0])u"hartree/bohr",
        (@SVector[0.0, 0.0, 0.0])u"hartree/bohr"
    ]
    @test virial_stress(system, p) isa SVector{6,<:Unitful.Energy}
    @test virial_stress(system, p) == (@SVector[0.0625, 0.0625, 0.0, 0.0, 0.0, 0.0])u"hartree"
end
