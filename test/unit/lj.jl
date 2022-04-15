@testset "Lennard-Jones Unit Tests" begin
    ϵ = 1.0u"eV"
    σ = 0.25u"bohr"
    rcutoff = 2.0u"Å"
    p = LennardJones(ϵ, σ, rcutoff, [:Ar, :H])

    @test p isa EmpiricalPotential
    @test p.ϵ == austrip(ϵ)
    @test p.σ == austrip(σ)
    @test p.rcutoff == austrip(rcutoff)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
