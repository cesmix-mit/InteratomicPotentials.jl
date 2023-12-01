push!(Base.LOAD_PATH, "../../")

using AtomsBase
using InteratomicPotentials
using LinearAlgebra
using Unitful
using UnitfulAtomic
using BenchmarkTools
using PProf

include("utils.jl")

@show N, NUM_ITER = parse.(Int, ARGS) # E.g. 100, 100

potential = LennardJones(1.0u"eV", 1.0u"Å", 100.0u"Å", [:Ar])

system = generate_lj_cluster(N)

display(@benchmark energy_and_force($system, $potential))

@pprof for _ = 1:NUM_ITER
    energy_and_force(system, potential)
end
