@testset "ZBL Unit Tests" begin
    Z₁ = 18
    Z₂ = 18
    e = 0.25u"C"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = ZBL(Z₁, Z₂, e, rcutoff, species)

    @test p isa EmpiricalPotential
    @test p.Z₁ == Z₁
    @test p.Z₂ == Z₂
    @test p.e == austrip(e)
    @test p.rcutoff == austrip(rcutoff)
    @test p.species == (:Ar, :H)

    @test get_parameters(p) == (;)
    @test set_parameters(p, (;)) === p
    @test serialize_parameters(p) == []
    @test deserialize_parameters(p, []) === p

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == ZBL(Z₁, Z₂, austrip(e), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == ZBL(Z₁, Z₂, austrip(e), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
