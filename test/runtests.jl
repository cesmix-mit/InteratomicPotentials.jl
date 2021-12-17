using Test

@time begin
@testset "Potentials.jl" begin
    # include("lj_test.jl")
    include("snap_test_single_element.jl")
end
end

