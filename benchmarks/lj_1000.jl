include("helpers.jl")

using BenchmarkTools

system = load_configuration("benchmarks/data/lj_clusters/1000.xyz")

@benchmark energy_and_force($system, $test_potential)
