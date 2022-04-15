@testset "Emperical Potential Default Implmentation Unit Tests" begin
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
    lj = LennardJones(ϵ * u"hartree", σ * u"bohr", rcutoff * u"bohr", [:Ar, :H])

    @test energy_and_force(system, lj).e isa ENERGY_TYPE
    @test energy_and_force(system, lj).f isa AbstractVector{SVector{3,FORCE_TYPE}}
    @test virial_stress(system, lj) isa SVector{6,ENERGY_TYPE}

    @test force(@SVector[1.0, 1.0, 1.0], lj) == force(sqrt(3), @SVector[1.0, 1.0, 1.0], lj)
end
