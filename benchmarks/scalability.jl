include("helpers.jl")

using BenchmarkTools
using Plots

ns = 250:250:10000
@time times = map(ns) do N
    system = generate_lj_cluster(N)
    time = @belapsed energy_and_force($system, $test_potential)
    @show N, time
    time
end

@show times

p = plot(title="Scaling of energy_and_force on LJ Clusters", xlabel="Number of Atoms", ylabel="Median Time Elapsed (s)", legend=false)
plot!(p, ns, times, markershape=:circle)
savefig(p, "scaling.svg")
