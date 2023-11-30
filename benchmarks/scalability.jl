push!(Base.LOAD_PATH, "../../")

using AtomsBase
using InteratomicPotentials
using LinearAlgebra
using Unitful
using UnitfulAtomic
using BenchmarkTools
using Plots

include("utils.jl")

potential = LennardJones(1.0u"eV", 1.0u"Å", 100.0u"Å", [:Ar])

ns = 250:250:1000
@time times = map(ns) do N
    system = generate_lj_cluster(N)
    time = @belapsed energy_and_force($system, $potential)
    @show N, time
    time
end

@show times

p = plot(title="Scaling of energy_and_force on LJ Clusters",
         xlabel="Number of Atoms",
         ylabel="Median Time Elapsed (s)",
         legend=false)
plot!(p, ns, times, markershape=:circle)
savefig(p, "scaling.svg")
