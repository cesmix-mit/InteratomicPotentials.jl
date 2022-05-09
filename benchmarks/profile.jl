include("helpers.jl")

using BenchmarkTools
using PProf

@show N, NUM_ITER = parse.(Int, ARGS)

potential = LennardJones(1.0u"eV", 1.0u"Å", 100.0u"Å", [:Ar])

system = generate_lj_cluster(N)

display(@benchmark energy_and_force($system, $potential))

@pprof for _ = 1:NUM_ITER
    energy_and_force(system, potential)
end
