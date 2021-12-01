struct SNA_BINDICES
    j1 :: Int
    j2 :: Int
    j  :: Int
end
struct SNA_ZINDICES 
    j1 :: Int
    j2 :: Int 
    j  :: Int
    ma1min :: Int 
    ma2max :: Int 
    mb1min :: Int 
    mb2max :: Int 
    na     :: Int 
    nb     :: Int 
    jju    :: Int
end
struct PrebuiltArrays{T}
    # Array indices, sizes
    idxcg_block :: Array{Int}
    idxu_block  :: Array{Int}
    idxb_block  :: Array{Int}
    idxz_block  :: Array{Int}
    idxcg_max :: Int
    idxu_max  :: Int
    idxb_max  :: Int
    idxz_max  :: Int

    idxb      :: Vector{SNA_BINDICES}
    idxz      :: Vector{SNA_ZINDICES}
    ncoeff    :: Int 


    # Clebsch Gordan
    rootpqarray :: Array{T}
    cglist      :: Vector{T}

    # Iterparticle Arrays
    rij         :: Array{T} 
    inside      :: Vector{T} 
    wj          :: Vector{T} 
    rcutij      :: Vector{T}
    element     :: Vector{Int} 
end

struct SNAPParams{S<:Integer, T<:AbstractFloat}
    n_atoms    :: S
    twojmax    :: S 
    n_elements :: S

    rcut       :: T
    rmin0      :: T
    rfac0      :: T

    chem_flag  :: Bool
    bzero_flag :: Bool 
    bnorm_flag :: Bool
    switch_flag :: Bool
    prebuilt_flag :: Bool 
    prebuilt_arrays :: PrebuiltArrays{T}
end

function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, prebuilt_flag::Bool, prebuilt_arrays::PrebuiltArrays{T} ) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams{S, T}(n_atoms, twojmax, n_elements, 
        rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, true, prebuilt_arrays)
    
end

function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, prebuilt_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
    if prebuilt_flag
        error("Cannot pass prebuilt_flag = true without also passing prebuilt arrays.")
    else
        return SNAPParams{S, T}(n_atoms, twojmax, n_elements, 
        rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, true, 
        initialize_prebuilt_arrays(twojmax, n_atoms, n_elements, chem_flag))
    end
end


function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
        rcut::T, rmin0::T, rfac0::T, 
        chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
        switch_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
        prebuilt_arrays = initialize_prebuilt_arrays()
        return SNAPParams(n_atoms, twojmax, n_elements, 
        rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, bzero_flag, false,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, false, false,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, false, false, false,
    true, false)
end


struct RuntimeArrays{T}

    # U Arrays: Real, Imaginary Parts
    ulisttot_r  :: Vector{T} 
    ulisttot_i  :: Vector{T} 
    dulist_r    :: Array{T} 
    dulist_i    :: Array{T}

    ulist_r_ij  :: Array{T} 
    ulist_i_ij  :: Array{T} 

    # Z Arrays: Real, Imaginary Parts
    zlist_r     :: Vector{T}
    zlist_i     :: Vector{T} 

    # Y Arrays: Real, Imaginary Parts
    ylist_r     :: Vector{T} 
    ylist_i     :: Vector{T} 

    # Bispectrum Arrays
    blist       :: Vector{T} 
    dblist      :: Array{T} 
end

# /* ---------------------------------------------------------------------- */
## Initialize

function initialize_prebuilt_arrays(twojmax::Int, n_atoms::Int, n_elements::Int, chem_flag::Bool)
    jdim = twojmax+1
    ncount = 0
    idxcg_count = 0
    idxu_count = 0
    idxb_count = 0
    idxz_count = 0

    idxcg_block = zeros(Int, jdim, jdim, jdim)
    idxu_block = zeros(Int, jdim)
    idxb_block = zeros(Int, jdim, jdim, jdim)
    idxz_block = zeros(Int, jdim, jdim, jdim)
    
    for j1 = 0:twojmax
        for j2 = 0:j1
            for j = (j1-j2):2:minimum([twojmax, j1+j2])
                if j >= j1 
                    ncount += 1

                    idxb_block[j1+1, j2+1, j+1] = idxb_count
                    idxb_count += 1
                    
                end

                #idxcg
                idxcg_block[j1+1, j2+1, j+1] = idxcg_count
                for m1 = 0:j1
                    for m2 = 0:j2 
                        idxcg_count += 1
                    end
                end

                #idxz
                mb = 0
                while 2*mb < j
                    for ma = 0:j;
                      idxz_count +=1
                    end
                    mb +=1;
                end
                
            end
        end

        #idxu_block 
        idxu_block[j1+1] = idxu_count
        for mb = 0:j1 
            for ma = 0:j1
                idxu_count += 1
            end
        end

    end
    
    idxu_max  = idxu_count
    idxcg_max = idxcg_count
    idxb_max = idxb_count
    idxz_max = idxz_count
    idxb = Vector{SNA_BINDICES}(undef, idxb_max)
    idxz = Vector{SNA_ZINDICES}(undef, idxz_max)

    ## Second loop for idxb
    idxb_count = 0
    idxz_count = 0
    for j1 = 0:twojmax
        for j2 = 0:j1
            j = j1 - j2
            while j <= minimum([twojmax, j1+j2])
                if j >= j1 
                    idxb_count += 1
                    idxb[idxb_count] = SNA_BINDICES(j1, j2, j)
                end

                idxz_block[j1+1, j2+1, j+1] = idxz_count;

                # // find right beta[jjb] entry
                # // multiply and divide by j+1 factors
                # // account for multiplicity of 1, 2, or 3
                mb = 0;
                while 2*mb < j
                    for ma = 0:j 
                        ma1min = maximum([0, (2 * ma - j - j2 + j1)/2])
                        ma2max = (2 * ma - j - (2*ma1min - j1) + j2) / 2
                        na     = minimum([j1, (2 * ma - j + j2 + j1) / 2]) - ma1min + 1 
                        mb1min = maximum([0, (2 * mb - j - j2 + j1) / 2])
                        mb2max = (2 * mb - j - (2*mb1min - j1) + j2)/2
                        nb     =  minimum([j1, (2*mb - j + j2 + j1) / 2]) - mb1min + 1

                        jju = idxu_block[j+1] + (j+1)*mb + ma 
                        idxz[idxz_count+1] = SNA_ZINDICES(j1, j2, j, 
                                            ma1min, ma2max, na, 
                                            mb1min, mb2max, nb, jju)
                        idxz_count += 1
                    end
                    mb += 1;
                end
                j+=2;
            end
        end
    end

    # Initialize Root array 
    rootpqarray = zeros(twojmax+2, twojmax+2)        # rootpqarray
    rootpqarray = init_rootpqarray(twojmax, rootpqarray)

    # Initialize Clebsch Gordan
    cglist      = zeros(idxcg_max)                  # cglist 
    cglist      = init_clebsch_gordan(twojmax, cglist)

    # Initialize ij Arrays
    rij         = zeros(n_atoms, 3)  # rij 
    inside      = zeros(n_atoms)     # inside 
    wj          = zeros(n_atoms)     # wj 
    rcutij      = zeros(n_atoms)     # rcutij 
    element     = ones(Int, n_atoms)     # element 

    
    # Calculate number of coefficients
    ndoubles = n_elements * n_elements
    ntriples = n_elements * n_elements * n_elements
    if chem_flag
        ncoeff = ncount * n_elements * n_elements * n_elements
    else
        ncoeff = ncount 
    end

    return PrebuiltArrays(
    idxcg_block,
    idxu_block,
    idxb_block,
    idxz_block,
    idxcg_max,
    idxu_max,
    idxb_max,
    idxz_max,
    idxb,
    idxz,
    ncoeff,
    rootpqarray,
    cglist,
    rij,
    inside,
    wj,
    rcutij,
    element)
end

function initialize_runtime_arrays(snap :: SNAPParams)
    n_doubles = snap.n_elements * snap.n_elements
    n_triples = n_doubles * snap.n_elements
    return RuntimeArrays(
        zeros(snap.prebuilt_arrays.idxu_max*snap.n_elements),   # ulisttot_r 
        zeros(snap.prebuilt_arrays.idxu_max*snap.n_elements),   # ulisttot_i 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_r 
        zeros(snap.prebuilt_arrays.idxu_max, 3),                # dulist_i 
        zeros(snap.n_atoms, snap.prebuilt_arrays.idxu_max),                       # ulist_r_ij 
        zeros(snap.n_atoms, snap.prebuilt_arrays.idxu_max),                       # ulist_i_ij
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_r 
        zeros(snap.prebuilt_arrays.idxz_max*n_doubles),         # zlist_i 
        zeros(snap.prebuilt_arrays.idxu_max*snap.n_elements),   # ylist_r 
        zeros(snap.prebuilt_arrays.idxu_max*snap.n_elements),   # ylist_i
        zeros(snap.prebuilt_arrays.idxb_max*n_triples),         # blist
        zeros(snap.prebuilt_arrays.idxb_max*n_triples, 3)       # dblist
    )
end

function deltacg(j1::Int, j2::Int, j::Int)
    sfaccg = factorial( Int( (j1+j2+j) / 2 + 1 ) )
    return sqrt( 
            factorial( Int( (j1+j2 -j)/2 ) ) * 
            factorial( Int( (j1-j2+j)/2 ) )  * 
            factorial( Int( (-j1+j2+j) / 2) ) / sfaccg
    )
end

function init_clebsch_gordan(twojmax::Int, cglist::Vector{T}) where T<: AbstractFloat
    idxcg_count = 0;
    for j1 = 0:twojmax
        for j2 = 0:j1
            for j = (j1-j2):2:minimum([twojmax, j1+j2])
                for m1 = 0:j1 
                    aa2 = 2 * m1 - j1;
    
                    for m2 = 0:j2   
                    bb2 = 2 * m2 - j2;
                    m = (aa2 + bb2 + j) / 2;
        
                        if (m < 0) || (m > j) 
                            cglist[idxcg_count+1] = 0.0;
                            idxcg_count+=1;
                            continue;
                        end
        
                        sum = 0.0;
                        min_z = maximum([0, maximum([-(j-j2+aa2)/2, -(j-j1-bb2)/2])])
                        max_z = minimum([ (j1 + j2 - j) / 2, minimum([ (j1-aa2)/2, (j2+bb2) / 2 ])])
                        for z = min_z:max_z
                            ifac = ( (z % 2)==0 ) ? -1 : 1;
                            sum += ifac /
                                (factorial(Int(z)) *
                                factorial( Int( (j1 + j2 - j) / 2 - z) ) *
                                factorial( Int( (j1 - aa2) / 2 - z) ) *
                                factorial( Int( (j2 + bb2) / 2 - z) ) *
                                factorial( Int( (j - j2 + aa2) / 2 + z) ) *
                                factorial( Int( (j - j1 - bb2) / 2 + z)) );
                        end
        
                        cc2 = 2 * m - j;
                        dcg = deltacg(j1, j2, j);
                        sfaccg = sqrt(factorial(Int( (j1 + aa2) / 2) )  *
                                    factorial( Int( (j1 - aa2) / 2) ) *
                                    factorial( Int( (j2 + bb2) / 2) ) *
                                    factorial( Int( (j2 - bb2) / 2) ) *
                                    factorial( Int( (j  + cc2) / 2) ) *
                                    factorial( Int( (j  - cc2) / 2) )*
                                    (j + 1));

                        cglist[idxcg_count+1] = sum * dcg * sfaccg;
                        idxcg_count+=1;
                    end
                end
            end
        end
    end
    return cglist
end

function init_rootpqarray(twojmax::Int, rootpqarray::Matrix{T}) where T<:AbstractFloat
    for p = 1:twojmax+1
        for q = 1:twojmax+1
            rootpqarray[p+1, q+1] = sqrt(p / q)
        end
    end
    return rootpqarray
end

function compute_sfac(snap::SNAPParams, r::AbstractFloat, rcut::AbstractFloat)
    if ~ snap.switch_flag
        sfac = 1.0 
    else
        if r <= snap.rmin0
            sfac = 1.0 
        elseif r > rcut
            sfac =  0.0
        else
            rcutfac = pi / (rcut - snap.rmin0)
            sfac = 0.5 * ( cos( (r - snap.rmin0) * rcutfac ) + 1.0 )
        end
    end
    return sfac
end

function compute_dsfac(snap::SNAPParams, r::AbstractFloat, rcut::AbstractFloat)
    if ~ snap.switch_flag
        dsfac = 1.0
    else
        if r <= snap.rmin0
            dsfac = 0.0 
        elseif rcut > rcut
            dsfac = 0.0
        else
            rcutfac = pi / (rcut - snap.rmin0)
            dsfac = -0.5 * ( sin( (r - snap.rmin0) * rcutfac ) * rcutfac )
        end
    end
    return dsfac
end
