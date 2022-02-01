using Distances
using NearestNeighbors
using StaticArrays

struct NeighborList
    j::Vector{Vector{Int}}
    R::Vector{Vector{Float64}}
    r::Vector{Vector{SVector{3,Float64}}}
end
Base.length(nn::NeighborList) = length(nn.j)

function get_distance(L::SVector{3,<:AbstractFloat}, x::SVector{3,<:AbstractFloat}, y::SVector{3,<:AbstractFloat})
    broadcast(L, x, y) do Li, xi, yi
        if Li == Inf
            xi - yi
        else
            d = mod(xi - yi, Li)
            2d < Li ? d : d - Li
        end
    end
end

function neighborlist(A::AbstractSystem{3}, rcutoff::Float64)
    # Convert Positions to Matrix for Ball tree
    X = [SVector{3}(ustrip.(p)) for p ∈ position(A)]

    # Create Metric for Periodic Boundary Conditions
    periodic = periodicity(A)
    L = @SVector [periodic[i] ? ustrip(bounding_box(A)[i][i]) : Inf for i = 1:3]
    d = Distances.PeriodicEuclidean(L)

    # Build Ball tree
    tree = BallTree(X, d)

    # Intialize empty vectors
    j = Vector{Vector{Int64}}(undef, length(X))              # Neighbors
    R = Vector{Vector{Float64}}(undef, length(X))            # Distances
    r = Vector{Vector{SVector{3,Float64}}}(undef, length(X)) # Positions

    # Fill vectors
    for n in 1:length(X)
        neighbors = filter(m -> m > n, inrange(tree, X[n], rcutoff, true))
        j[n] = zeros(Int64, length(neighbors))
        R[n] = zeros(Float64, length(neighbors))
        r[n] = zeros(SVector{3,Float64}, length(neighbors))
        for (i, m) in enumerate(neighbors)
            rr = get_distance(L, X[n], X[m])
            j[n][i] = m
            R[n][i] = norm(rr)
            r[n][i] = rr
        end
    end
    NeighborList(j, R, r)
end
