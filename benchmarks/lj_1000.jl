include("helpers.jl")

using BenchmarkTools

T, = parse.(Int, ARGS)

potential = LennardJones(1.0u"eV", 1.0u"Å", 100.0u"Å", [:Ar])

system = load_configuration("benchmarks/data/lj_clusters/1000.xyz")

display(@benchmark energy_and_force($system, $potential) seconds = T)
