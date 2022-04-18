@testset "Born-Mayer Unit Tests" begin
    A = 1.0u"eV"
    ρ = 0.25u"bohr"
    σ = 0.25u"bohr"
    C = 1.0u"eV*Å"
    D = 1.0u"eV*Å"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = BornMayer(A, ρ, σ, C, D, rcutoff, species)

    @test p isa EmpiricalPotential{NamedTuple{(:A, :ρ, :σ, :C, :D)},NamedTuple{(:rcutoff,)}}

    @test get_rcutoff(p) == austrip(rcutoff)
    @test get_species(p) == (:Ar, :H)

    @test get_parameters(p) == (; A=austrip(A), ρ=austrip(ρ), σ = austrip(σ), C = austrip(C), D = austrip(D))
    @test set_parameters(p, (; A=2.0, ρ=3.0, σ = 4.0, C = 5.0, D = 6.0)) == BornMayer(2.0, 3.0, 4.0, 5.0, 6.0, austrip(rcutoff), species)
    @test serialize_parameters(p) == [austrip(A), austrip(ρ), austrip(σ), austrip(C), austrip(D)]
    @test deserialize_parameters(p, [2.0, 3.0, 4.0, 5.0, 6.0]) == BornMayer(2.0, 3.0, 4.0, 5.0, 6.0, austrip(rcutoff), species)

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == BornMayer(austrip(A), austrip(ρ), austrip(σ), austrip(C), austrip(D), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == BornMayer(austrip(A), austrip(ρ), austrip(σ), austrip(C), austrip(D), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
