push!(Base.LOAD_PATH, "../../")

using AtomsBase
using InteratomicPotentials
using LinearAlgebra
using Unitful
using UnitfulAtomic
using BenchmarkTools
using Plots
using BenchmarkTools

include("utils.jl")

T, = parse.(Int, ARGS)

potential = LennardJones(1.0u"eV", 1.0u"Å", 100.0u"Å", [:Ar])

system = load_configuration("data/lj_clusters/1000.xyz")

display(@benchmark energy_and_force($system, $potential) seconds = T)
