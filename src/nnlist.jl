# ! The contents of this file will eventually be made obsolete by the creation of a general neighbor list package.

export NeighborList, neighborlist

struct NeighborList
    j::Vector{Vector{Int}}                 # Neighbors
    r::Vector{Vector{SVector{3,Float64}}}  # Displacements
    function NeighborList(n)
        new(
            Vector{Vector{Int64}}(undef, n),
            Vector{Vector{SVector{3,Float64}}}(undef, n)
        )
    end
end
Base.length(nn::NeighborList) = length(nn.j)

# calculates displacement vector x - y with respect to the boundary conditions L
function get_displacement(L::SVector{3,F}, x::SVector{3,F}, y::SVector{3,F}) where {F<:AbstractFloat}
    d = mod.(x - y, L)
    @. ifelse(isinf(L), x - y, ifelse(2d < L, d, d - L))
end

function neighborlist(s::AbstractSystem{3}, rcutoff::Float64)
    # Convert Positions to Matrix for Ball tree
    X = [SVector{3}(austrip.(p)) for p ∈ position(s)]::Vector{SVector{3,Float64}}

    # Create Metric for Periodic Boundary Conditions
    periodic = periodicity(s)
    L = @SVector [periodic[i] ? austrip(bounding_box(s)[i][i]) : Inf for i ∈ 1:3]
    d = Distances.PeriodicEuclidean(L)

    # Build Ball tree
    tree = BallTree(X, d)

    # Intialize empty vectors
    nl = NeighborList(length(X))

    # Fill vectors
    @inbounds for n in 1:length(X)
        nl.j[n] = filter(m -> m > n, inrange(tree, X[n], rcutoff))
        nl.r[n] = map(nl.j[n]) do m
            get_displacement(L, X[n], X[m])
        end
    end
    nl
end
