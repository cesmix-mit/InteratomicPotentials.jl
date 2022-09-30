using Test

@testset "InteratomicBasisPotentials.jl" begin

    @testset "SNAP" begin
        include("SNAP/snap_test_single_element.jl")
        include("SNAP/snap_test_multi_element.jl")
        # include("SNAP/snap_test_performance.jl")
    end

    @testset "ACE" begin
        include("ACE/ace_test.jl")
    end
end
