@testset "Coulomb Unit Tests" begin
    q₁ = 1.0u"C"
    q₂ = 1.0u"C"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = Coulomb(q₁, q₂, rcutoff, species)

    @test p isa EmpiricalPotential
    @test p.q₁ == austrip(q₁)
    @test p.q₂ == austrip(q₂)
    @test p.rcutoff == austrip(rcutoff)
    @test p.species == (:Ar, :H)

    @test get_parameters(p) == (;)
    @test set_parameters(p, (;)) === p
    @test serialize_parameters(p) == []
    @test deserialize_parameters(p, []) === p

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == Coulomb(austrip(q₁), austrip(q₂), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == Coulomb(austrip(q₁), austrip(q₂), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
