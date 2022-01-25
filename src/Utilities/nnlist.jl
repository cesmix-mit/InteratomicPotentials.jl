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

to_array(A::AbstractSystem) = [SVector{3}(ustrip.(p)) for p ∈ position(A)]

function neighborlist(A::AbstractSystem{3}, rcutoff::Float64)
    # Convert Positions to Matrix for Ball tree
    X = to_array(A)

    # Create Metric for Periodic Boundary Conditions
    periodic = periodicity(A)
    L = @SVector [periodic[i] ? ustrip(bounding_box(A)[i][i]) : Inf for i = 1:3]
    d = Distances.PeriodicEuclidean(L)

    # Build Ball tree
    tree = BallTree(X, d)

    # Intialize empty vectors
    j = Vector{Int64}[] # j 
    R = Vector{Float64}[] # Distances
    r = Vector{SVector{3,Float64}}[] # Positions

    # Fill vectors
    for n = 1:length(X)
        neighbors = inrange(tree, X[n], rcutoff, true)
        jtemp = Int64[]
        Rtemp = Float64[]
        rtemp = SVector{3,Float64}[]
        for neighbor in neighbors
            if neighbor != n
                push!(jtemp, neighbor)
                rr = get_distance(L, X[n], X[neighbor])
                push!(rtemp, rr)
                push!(Rtemp, norm(rr))
            end
        end
        push!(j, jtemp)
        push!(R, Rtemp)
        push!(r, rtemp)
    end
    NeighborList(j, R, r)
end
