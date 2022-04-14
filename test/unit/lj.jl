@testset "Lennard-Jones" begin
    position = (@SVector [1.0, 0.0, 0.0])u"bohr"
    atoms = [
        :Ar => 0.5 * position,
        :Ar => 0.75 * position,
        :H => 1.1 * position
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    ϵ = 1.0
    σ = 0.25
    rcutoff = 2.0
    lj = LennardJones(ϵ * u"hartree", σ * u"bohr", rcutoff * u"bohr", [:Ar, :H])

    @test lj isa EmpiricalPotential
    @test potential_energy(system, lj) isa ENERGY_TYPE
    @test force(system, lj) isa AbstractVector{SVector{3,FORCE_TYPE}}
    @test virial(system, lj) isa ENERGY_TYPE
    @test virial_stress(system, lj) isa SVector{6,ENERGY_TYPE}
end