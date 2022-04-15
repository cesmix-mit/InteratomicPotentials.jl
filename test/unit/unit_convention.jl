import Unitful.𝐌, Unitful.𝐋, Unitful.𝐓

@testset "Unit Convention" begin
    @test InteratomicPotentials.ENERGY_UNIT == aunit(𝐌 * 𝐋^2 / 𝐓^2)
    @test InteratomicPotentials.FORCE_UNIT == aunit(𝐌 * 𝐋 / 𝐓^2)

    @test InteratomicPotentials.ENERGY_TYPE == typeof(1.0aunit(𝐌 * 𝐋^2 / 𝐓^2))
    @test InteratomicPotentials.FORCE_TYPE == typeof(1.0aunit(𝐌 * 𝐋 / 𝐓^2))
end
