@testset "ZBL Unit Tests" begin
    Z₁ = 18
    Z₂ = 18
    e = 0.25u"C"
    rcutoff = 2.0u"Å"
    p = ZBL(Z₁, Z₂, e, rcutoff, [:Ar, :H])

    @test p isa EmpiricalPotential
    @test p.Z₁ == austrip(Z₁)
    @test p.Z₂ == austrip(Z₂)
    @test p.e == austrip(e)
    @test p.rcutoff == austrip(rcutoff)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
