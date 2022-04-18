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
    @testset "Unit Tests" begin
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
    @testset "Integration Tests" begin
        # TODO: Dallas
    end
end
