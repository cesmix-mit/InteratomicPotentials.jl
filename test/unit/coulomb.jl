@testset "Coulomb Unit Tests" begin
    q₁ = 1.0u"C"
    q₂ = 1.0u"C"
    rcutoff = 2.0u"Å"
    p = Coulomb(q₁, q₂, rcutoff, [:Ar, :H])

    @test p isa EmpiricalPotential
    @test p.q₁ == austrip(q₁)
    @test p.q₂ == austrip(q₂)
    @test p.rcutoff == austrip(rcutoff)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
