using Test

@time begin
@testset "Potentials.jl" begin
    include("lj_test.jl")
end
end
