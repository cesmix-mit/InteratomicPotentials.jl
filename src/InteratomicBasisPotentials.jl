module InteratomicBasisPotentials

using InteratomicPotentials
using StaticArrays
using LinearAlgebra
using AtomsBase
using Unitful, UnitfulAtomic
import Flux: gradient, Chain, Dense
import Zygote: withgradient

include("api.jl")
include("types.jl")


end # module
