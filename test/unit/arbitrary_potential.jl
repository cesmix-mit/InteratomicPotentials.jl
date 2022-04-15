@testset "Arbitrary Potential Default Implmentation Unit Tests" begin
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

    @test potential_energy(system, lj) isa ENERGY_TYPE
    @test potential_energy(system, lj) == energy_and_force(system, lj).e
    @test force(system, lj) isa AbstractVector{SVector{3,FORCE_TYPE}}
    @test force(system, lj) == energy_and_force(system, lj).f
    @test virial(system, lj) isa ENERGY_TYPE
    @test virial(system, lj) == sum(virial_stress(system, lj)[1:3])
end
