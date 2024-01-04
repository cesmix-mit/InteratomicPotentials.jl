using InteratomicPotentials

using AtomsBase
using LinearAlgebra
using StaticArrays
using Test
using Unitful
using UnitfulAtomic

include("mocks.jl")

@testset "Potentials.jl" begin
    # The first set of tests are simple unit tests that generally test just one module at a time.
    # These tests collectively have high code coverage and test whether each module is structurally correct.
    # They do not validate the numeric results for scientific accuracy.
    @time @testset "Unit Tests" begin
        include("unit/constants.jl")
        include("unit/unit_convention.jl")
        include("unit/nnlist.jl")
        @testset "Default Implementations" begin
            include("unit/types/abstract_potential.jl")
            include("unit/types/non_trainable_potential.jl")
            include("unit/types/trainable_potential.jl")
            include("unit/types/empirical_potential.jl")
        end
        @testset "Empirical Potentials" begin
            include("unit/empirical/lj.jl")
            include("unit/empirical/bm.jl")
            include("unit/empirical/coulomb.jl")
            include("unit/empirical/zbl.jl")
            include("unit/empirical/morse.jl")
        end
        @testset "Mixed Potentials" begin
            include("unit/mixed/linear_combination_potential.jl")
        end
    end

    # The second set of tests are integration tests which incorporate multiple modules.
    # These tests validate the numeric results of calcuations to ensure that the modules are logically correct.
    # In other words, they test that the full system works together and produces scientifically accurate results.
    @time @testset "Integration Tests" begin
        include("integration/lj_clusters/lj_150.jl")
        include("integration/lj_clusters/lj_1000.jl")
    end

    @time @testset "SNAP" begin
        include("SNAP/snap_test_single_element.jl")
        include("SNAP/snap_test_multi_element.jl")
        # include("SNAP/snap_test_performance.jl")
    end

    @time @testset "ACE" begin
        include("ACE/ace_test.jl")
    end
end
