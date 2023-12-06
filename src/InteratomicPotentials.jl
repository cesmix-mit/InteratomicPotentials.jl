module InteratomicPotentials
using Requires  # Only needed in Julia < 1.9

using AtomsBase
using Base.Threads
using LinearAlgebra
using StaticArrays
using Unitful
using UnitfulAtomic

using Distances
using NearestNeighbors

include("unit_convention.jl")
include("constants.jl")
include("nnlist.jl")
include("api.jl")
include("types.jl")

# Requires-based dependency management for Julia < 1.9
if !isdefined(Base, :get_extension)
    function __init__()
        @require DFTK="acf6eb54-70d9-11e9-0013-234b7a5f5337" begin
            include("../ext/InteratomicPotentialsDFTKExt.jl")
        end
    end
end

end
