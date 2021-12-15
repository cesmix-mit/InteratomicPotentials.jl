struct SNAPParams{S<:Integer, T<:AbstractFloat}
    n_atoms    :: S
    twojmax    :: S 
    elements   :: Vector{Symbol}

    rcut       :: Vector{T}
    rmin0      :: T
    rfac0      :: T
    weight     :: Vector{T}

    bzero_flag :: Bool 
    bnorm_flag :: Bool
    switch_flag :: Bool
    wselfall_flag :: Bool
    prebuilt_flag :: Bool 
    prebuilt_arrays :: PrebuiltArrays{T}
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{T}, rmin0::T, rfac0::T, weight ::Vector{T},
    bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, prebuilt_flag::Bool, prebuilt_arrays::PrebuiltArrays{T} ) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcut, rmin0, rfac0, weight, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, true, prebuilt_arrays)
    
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T},
    bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, wselfall_flag::Bool, prebuilt_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
    if prebuilt_flag
        error("Cannot pass prebuilt_flag = true without also passing prebuilt arrays.")
    else
        return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcut, rmin0, rfac0, weight, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, true, 
        initialize_prebuilt_arrays(twojmax, length(elements), bzero_flag, bnorm_flag))
    end
end


function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
        rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T},
        bzero_flag::Bool, bnorm_flag::Bool,
        switch_flag::Bool, wselfall_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
        return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcut, rmin0, rfac0, weight, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, false)
end
function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T},
    bzero_flag::Bool, bnorm_flag::Bool, switch_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcut, rmin0, rfac0, weight, bzero_flag, bnorm_flag,
    true, false, false)
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T},
    bzero_flag::Bool, bnorm_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcut, rmin0, rfac0, weight, bzero_flag, bnorm_flag,
    true, false, false)
end
function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T},
    bzero_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcut, rmin0, rfac0, weight, bzero_flag, false,
    true, false, false)
end
function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcut::Vector{<:T}, rmin0::T, rfac0::T,  weight ::Vector{<:T}) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcut, rmin0, rfac0, weight, false, false,
    false, false, false)
end
