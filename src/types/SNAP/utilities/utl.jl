using Distances
using NearestNeighbors
include("compute_ui.jl")
include("compute_z.jl")
include("compute_b.jl")
include("compute_yi.jl")
include("compute_dij.jl")
function get_num_snap_coeffs(twojmax :: Int, n_elems :: Int, chem_flag :: Bool)
    J = twojmax 
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer!")
    end 
    if chem_flag
        num_coeffs = num_coeffs * n_elems * n_elems * n_elems
    else
        num_coeffs = num_coeffs 
    end
    return num_coeffs 
end

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

function neighborlist(A::AbstractSystem, snap::SNAPParams)

    # Construct cutoff
    num_elements = length(snap.species) 
    cutsq = zeros(num_elements, num_elements)
    cutmax = 0.0
    for i = 1:num_elements
        cut = 2.0*snap.radii[i]*snap.rcutfac
        if cut > cutmax
            cutmax = cut
        end
        cutsq[i, i] = cut*cut

        for j = (i+1):num_elements
            cut = (snap.radii[i] + snap.radii[j]) * snap.rcutfac
            cutsq[i, j] = cut*cut
            cutsq[j, i] = cut*cut
        end
    end
        

    # Convert Positions to Matrix for Ball tree
    X = [SVector{3}(ustrip.(p)) for p ∈ position(A)]

    # Create Metric for Periodic Boundary Conditions
    periodic = periodicity(A)
    L = @SVector [ periodic[i] ? ustrip(bounding_box(A)[i][i]) : Inf for i = 1:3]
    d = Distances.PeriodicEuclidean(L)
    
    # Build Ball tree
    tree = BallTree(X, d)

    # elements = unique([element(ai) for ai in A])
    elements = unique(atomic_symbol(A))

    # Intialize empty vectors
    j = Vector{Int64}[] # j 
    R = Vector{Float64}[] # Distances
    r = Vector{SVector{3, Float64}}[] # Positions

    # Fill vectors
    for n = 1:length(X)
        n_element = findall(x->x==atomic_symbol(A, n), elements)[1]
        neighbors = inrange(tree, X[n], sqrt(cutmax), true)
        jtemp = Int64[]
        Rtemp = Float64[]
        rtemp = SVector{3, Float64}[]
        for (m, neighbor) in enumerate(neighbors)
            m_element = findall(x->x==atomic_symbol(A, m), elements)[1]
            rr = -1.0get_distance(L, X[n], X[neighbor])
            rsq = rr[1]*rr[1] + rr[2]*rr[2] + rr[3]*rr[3]
            if (neighbor != n) & (rsq <= cutsq[n_element, m_element])
                push!(jtemp, neighbor)
                push!(rtemp, rr)
                push!(Rtemp, sqrt(rsq))
            end
        end
    
        push!(j, jtemp)
        push!(R, Rtemp)
        push!(r, rtemp)
    end
    return NeighborList(j, R, r)
end
