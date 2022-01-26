using Test

@time @testset "Potentials.jl" begin
    include("lj_test.jl")
    include("SNAP/snap_test_single_element.jl")
    include("SNAP/snap_test_multi_element.jl")
    include("ace_test.jl")
end
