@testset "Born-Mayer Unit Tests" begin
    A = 1.0u"eV"
    ρ = 0.25u"bohr"
    rcutoff = 2.0u"Å"
    p = BornMayer(A, ρ, rcutoff, [:Ar, :H])

    @test p isa EmpiricalPotential
    @test p.A == austrip(A)
    @test p.ρ == austrip(ρ)
    @test p.rcutoff == austrip(rcutoff)

    r = @SVector[1.0, 1.0, 1.0]
    R = norm(r)

    @test potential_energy(R, p) isa Float64
    @test force(R, r, p) isa SVector{3,Float64}
end
