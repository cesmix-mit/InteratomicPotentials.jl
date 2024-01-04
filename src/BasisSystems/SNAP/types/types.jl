## TODO: Documentation
struct NeighborListSNAP 
    j 
    R 
    r 
end
include("prebuilt_arrays.jl")
struct SNAP{S<:Integer, T<:AbstractFloat} <: BasisSystem
    n_atoms    :: S
    twojmax    :: S 
    species   :: Vector{Symbol}

    rcutfac    :: T
    rmin0      :: T
    rfac0      :: T
    radii      :: Vector{T}
    weight     :: Vector{T}
    chem_flag  :: Bool
    bzero_flag :: Bool 
    bnorm_flag :: Bool
    switch_flag :: Bool
    wselfall_flag :: Bool
    prebuilt_flag :: Bool 
    prebuilt_arrays :: PrebuiltArrays{T}
end
include("runtime_arrays.jl")