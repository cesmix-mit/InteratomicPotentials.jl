using Test

@time @testset "Potentials.jl" begin
    include("lj_test.jl")
    include("suttonchen.jl")
end
