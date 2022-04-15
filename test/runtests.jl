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
            include("unit/arbitrary_potential.jl")
            include("unit/empirical_potential.jl")
        end
        @testset "Empirical Potentials" begin
            include("unit/lj.jl")
            include("unit/bm.jl")
            include("unit/coulomb.jl")
            include("unit/zbl.jl")
            include("unit/morse.jl")
        end
        @testset "Mixed Potentials" begin
            include("unit/mixed_potentials.jl")
        end
    end

    # The second set of tests are integration tests which incorporate multiple modules.
    # These tests validate the numeric results of calcuations to ensure that the modules are logically (scientifically) correct.
    @testset "Integration Tests" begin
        # TODO: Dallas
    end
end
