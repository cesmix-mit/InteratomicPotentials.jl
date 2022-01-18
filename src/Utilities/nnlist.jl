using Distances
using NearestNeighbors
using StaticArrays

struct NeighborList
    i::Vector{Vector{Int}}
    j::Vector{Vector{Int}}
    R::Vector{Vector{Float64}}
    r::Vector{Vector{SVector{3,Float64}}}
end
function length(nn::NeighborList)
    return length(nn.i)
end

function get_distance(L::SVector{3,<:AbstractFloat}, x::SVector{3,<:AbstractFloat}, y::SVector{3,<:AbstractFloat})
    r = [0.0, 0.0, 0.0]
    for i = 1:3
        if L == Inf
            r[i] += y[i] - x[i]
        else
            m = mod(abs(y[i] - x[i]), L[i])
            if m < L[i] - m
                r[i] += y[i] - x[i]
            else
                r[i] += L[i] - (y[i] - x[i])
            end
        end
    end
    return SVector{3}(r)
end

function to_array(A::AbstractSystem)
    positions = position(A)
    X = SVector{3,Float64}[]
    for p in positions
        push!(X, ustrip.(p))
    end
    return X
end
function neighborlist(A::AbstractSystem, rcutoff::Float64)

    # Convert Positions to Matrix for Ball tree
    X = to_array(A)

    # Create Metric for Periodic Boundary Conditions
    periodic = periodicity(A)
    L = @SVector [periodic[i] ? ustrip(bounding_box(A)[i][i]) : Inf for i = 1:3]
    d = Distances.PeriodicEuclidean(L)

    # Build Ball tree
    tree = BallTree(X, d)

    # Intialize empty vectors
    i = Vector{Int64}[] # i 
    j = Vector{Int64}[] # j 
    R = Vector{Float64}[] # Distances
    r = Vector{SVector{3,Float64}}[] # Positions

    # Fill vectors
    for n = 1:length(X)
        neighbors = inrange(tree, X[n], rcutoff, true)
        itemp = Int64[]
        jtemp = Int64[]
        Rtemp = Float64[]
        rtemp = SVector{3,Float64}[]
        for m in 1:length(neighbors)
            if neighbors[m] != n
                push!(itemp, n)
                push!(jtemp, neighbors[m])
                rr = get_distance(L, X[n], X[neighbors[m]])
                push!(rtemp, rr)
                push!(Rtemp, norm(rr))
            end
        end

        push!(i, itemp)
        push!(j, jtemp)
        push!(R, Rtemp)
        push!(r, rtemp)
    end
    return NeighborList(i, j, R, r)
end
