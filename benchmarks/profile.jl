include("helpers.jl")

using BenchmarkTools
using PProf

@show N, NUM_ITER = parse.(Int, ARGS)

system = generate_lj_cluster(N)

display(@benchmark energy_and_force($system, $test_potential))

@pprof for _ = 1:NUM_ITER
    energy_and_force(system, test_potential)
end
