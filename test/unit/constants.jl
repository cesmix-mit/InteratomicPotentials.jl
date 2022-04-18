@testset "Physical Constants" begin
    @test InteratomicPotentials.ϵ₀ ≈ austrip(8.8541878128e-12u"F / m")
    @test InteratomicPotentials.kₑ ≈ austrip(8.9875517923e9u"N * m^2 / C^2")
end
