struct RuntimeArrays{T}
    # Iterparticle Arrays
    indij       :: Vector{SVector{2, Int}}
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

function initialize_runtime_arrays(i::Int, nnlist::NeighborList, A::AbstractSystem, snap :: SNAPParams)
    # Initialize Arrays
    ind_ij  = SVector{2, Int}[]  # Contains (i,j) pairs 
    rcut_ij = AbstractFloat[]   # Contains rcut between (i,j)
    inside  = AbstractFloat[]   # Unused 
    wj      = AbstractFloat[]   # Weight of particle (j) w.r.t particle (i)
    rij     = SVector{3, AbstractFloat}[] # vector from particle j to particle i 
    element = Int[]            # Element of particle (j) 
    elements = snap.elements

    i_list = nnlist.i[i]
    j = nnlist.j[i]
    r = nnlist.r[i]
    i_element = findall(x->x==Symbol(A[i].atomic_symbol), elements)[1]
    for (jj, rj) in zip(j, r)
        push!(ind_ij, SVector{2}([i, jj]))

        j_element = findall(x->x==Symbol(A[jj].atomic_symbol), elements)[1]
        push!(element, j_element)
        push!(rcut_ij, (snap.radii[i_element] + snap.radii[j_element])*snap.rcutfac)
        push!(wj, snap.weight[j_element])
        push!(inside, 1)
        push!(rij, rj)
    end


    n_interactions = length(rij)
    num_elements = snap.prebuilt_arrays.num_elements
    n_doubles = num_elements * num_elements
    n_triples = n_doubles * num_elements
    return RuntimeArrays{AbstractFloat}(
        ind_ij,                                                 # indij
        rij,                                                    # rij
        inside,                                                 # inside
        wj,                                                     # wj
        rcut_ij,                                                 # rcutij
        element,                                                # element
        zeros(snap.prebuilt_arrays.idxu_max*num_elements),   # ulisttot_r 
        zeros(snap.prebuilt_arrays.idxu_max*num_elements),   # ulisttot_i 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_r 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_i 
        zeros(n_interactions+1, snap.prebuilt_arrays.idxu_max),                       # ulist_r_ij 
        zeros(n_interactions+1, snap.prebuilt_arrays.idxu_max),                       # ulist_i_ij
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_r 
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_i 
        zeros(snap.prebuilt_arrays.idxu_max*num_elements),   # ylist_r 
        zeros(snap.prebuilt_arrays.idxu_max*num_elements),   # ylist_i
        zeros(snap.prebuilt_arrays.idxb_max*n_triples),         # blist
        zeros(snap.prebuilt_arrays.idxb_max*n_triples, 3)       # dblist
    )
end