@testset "Born-Mayer Unit Tests" begin
    A = 1.0u"eV"
    ρ = 0.25u"bohr"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = BornMayer(A, ρ, rcutoff, species)

    @test p isa EmpiricalPotential
    @test p.A == austrip(A)
    @test p.ρ == austrip(ρ)
    @test p.rcutoff == austrip(rcutoff)
    @test p.species == (:Ar, :H)

    @test get_parameters(p) == (; A=austrip(A), ρ=austrip(ρ))
    @test set_parameters(p, (; A=2.0, ρ=3.0)) == BornMayer(2.0, 3.0, austrip(rcutoff), species)
    @test serialize_parameters(p) == [austrip(A), austrip(ρ)]
    @test deserialize_parameters(p, [2.0, 3.0]) == BornMayer(2.0, 3.0, austrip(rcutoff), species)

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == BornMayer(austrip(A), austrip(ρ), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == BornMayer(austrip(A), austrip(ρ), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
