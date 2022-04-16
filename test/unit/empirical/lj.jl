@testset "Lennard-Jones Unit Tests" begin
    ϵ = 1.0u"eV"
    σ = 0.25u"bohr"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = LennardJones(ϵ, σ, rcutoff, species)

    @test p isa EmpiricalPotential
    @test p.ϵ == austrip(ϵ)
    @test p.σ == austrip(σ)
    @test p.rcutoff == austrip(rcutoff)
    @test p.species == (:Ar, :H)

    @test get_parameters(p) == (; ϵ=austrip(ϵ), σ=austrip(σ))
    @test set_parameters(p, (; ϵ=2.0, σ=3.0)) == LennardJones(2.0, 3.0, austrip(rcutoff), species)
    @test serialize_parameters(p) == [austrip(ϵ), austrip(σ)]
    @test deserialize_parameters(p, [2.0, 3.0]) == LennardJones(2.0, 3.0, austrip(rcutoff), species)

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == LennardJones(austrip(ϵ), austrip(σ), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == LennardJones(austrip(ϵ), austrip(σ), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
