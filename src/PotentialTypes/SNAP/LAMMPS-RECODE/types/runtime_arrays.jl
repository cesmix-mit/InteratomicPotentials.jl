
struct RuntimeArrays{T}

    # Iterparticle Arrays
    rij         :: Vector{SVector{3, <:T}} 
    inside      :: Vector{<:T} 
    wj          :: Vector{<:T} 
    rcutij      :: Vector{<:T}
    element     :: Vector{Int} 

    # U Arrays: Real, Imaginary Parts
    ulisttot_r  :: Vector{<:T} 
    ulisttot_i  :: Vector{<:T} 
    dulist_r    :: Matrix{<:T} 
    dulist_i    :: Matrix{<:T}

    ulist_r_ij  :: Matrix{<:T} 
    ulist_i_ij  :: Matrix{<:T} 

    # Z Arrays: Real, Imaginary Parts
    zlist_r     :: Vector{<:T}
    zlist_i     :: Vector{<:T} 

    # Y Arrays: Real, Imaginary Parts
    ylist_r     :: Vector{<:T} 
    ylist_i     :: Vector{<:T} 

    # Bispectrum Arrays
    blist       :: Vector{<:T} 
    dblist      :: Matrix{<:T} 
end


function initialize_runtime_arrays(a::StaticAtom, A::AbstractSystem, snap :: SNAPParams)
    rcutij     = AbstractFloat[] 
    inside  = AbstractFloat[]
    wj      = AbstractFloat[]
    rij  = SVector{3, AbstractFloat}[] 
    element = Int[]
    elements = snap.elements
    for p_j = A.particles
        elem_temp = findall(x->x==Symbol(p_j.element.symbol), elements)[1]
        rcutij_temp = snap.rcut[elem_temp]
        wj_temp     = snap.weight[elem_temp]
        rij_temp = ustrip.(a.position - p_j.position )
        norm_ij = norm(rij_temp)
        if (norm_ij < snap.rmin0) || (norm_ij > rcutij_temp)
            continue;
        else
            push!(element, elem_temp)
            push!(rcutij, rcutij_temp)
            push!(wj, wj_temp)
            push!(inside, 1)
            push!(rij, rij_temp)
        end
    end
    

    n_interactions = length(rij)
    n_elements = length(snap.elements)
    n_doubles = n_elements * n_elements
    n_triples = n_doubles * n_elements
    return RuntimeArrays{AbstractFloat}(
        rij,                                                    # rij
        inside,                                                 # inside
        wj,                                                     # wj
        rcutij,                                                 # rcutij
        element,                                                # element
        zeros(snap.prebuilt_arrays.idxu_max*n_elements),   # ulisttot_r 
        zeros(snap.prebuilt_arrays.idxu_max*n_elements),   # ulisttot_i 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_r 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_i 
        zeros(n_interactions, snap.prebuilt_arrays.idxu_max),                       # ulist_r_ij 
        zeros(n_interactions, snap.prebuilt_arrays.idxu_max),                       # ulist_i_ij
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_r 
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_i 
        zeros(snap.prebuilt_arrays.idxu_max*n_elements),   # ylist_r 
        zeros(snap.prebuilt_arrays.idxu_max*n_elements),   # ylist_i
        zeros(snap.prebuilt_arrays.idxb_max*n_triples),         # blist
        zeros(snap.prebuilt_arrays.idxb_max*n_triples, 3)       # dblist
    )
end