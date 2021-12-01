struct SNAPParams{S, T} where S<:Integer, T<:AbstractFloat
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
        switch_flag::Bool ) where S<:Integer, T<:AbstractFloat
        return SNAPParams(n_atoms, twojmax, n_elements, 
        rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool) where S<:Integer, T<:AbstractFloat
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, bzero_flag, bnorm_flag,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool, bzero_flag::Bool) where S<:Integer, T<:AbstractFloat
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, bzero_flag, false,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T, 
    chem_flag::Bool) where S<:Integer, T<:AbstractFloat
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, chem_flag, false, false,
    true, false)
end
function SNAPParams(n_atoms::S, twojmax::S, n_elements::S, 
    rcut::T, rmin0::T, rfac0::T) where S<:Integer, T<:AbstractFloat
    return SNAPParams(n_atoms, twojmax, n_elements, 
    rcut, rmin0, rfac0, false, false, false,
    true, false)
end


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
    element     :: Vector{T} 
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