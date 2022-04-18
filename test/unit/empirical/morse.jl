@testset "Morse Unit Tests" begin
    D = 1.0u"eV"
    α = 2.0
    σ = 0.25u"bohr"
    rcutoff = 2.0u"Å"
    species = [:Ar, :H]
    p = Morse(D, α, σ, rcutoff, species)

    @test p isa EmpiricalPotential{NamedTuple{(:D, :α, :σ)},NamedTuple{(:rcutoff,)}}

    @test get_rcutoff(p) == austrip(rcutoff)
    @test get_species(p) == (:Ar, :H)

    @test get_parameters(p) == (; D=austrip(D), α, σ=austrip(σ))
    @test set_parameters(p, (; D=2.0, α=1.0, σ=3.0)) == Morse(2.0, 1.0, 3.0, austrip(rcutoff), species)
    @test serialize_parameters(p) == [austrip(D), α, austrip(σ)]
    @test deserialize_parameters(p, [2.0, 1.0, 3.0]) == Morse(2.0, 1.0, 3.0, austrip(rcutoff), species)

    @test get_hyperparameters(p) == (; rcutoff=austrip(rcutoff))
    @test set_hyperparameters(p, (; rcutoff=1.0)) == Morse(austrip(D), α, austrip(σ), 1.0, species)
    @test serialize_hyperparameters(p) == [austrip(rcutoff)]
    @test deserialize_hyperparameters(p, [1.0]) == Morse(austrip(D), α, austrip(σ), 1.0, species)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
