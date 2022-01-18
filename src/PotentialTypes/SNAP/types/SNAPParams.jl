struct SNAPParams{S<:Integer, T<:AbstractFloat}
    n_atoms    :: S
    twojmax    :: S 
    elements   :: Vector{Symbol}

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

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, prebuilt_flag::Bool, prebuilt_arrays::PrebuiltArrays{T} ) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, true, prebuilt_arrays)
    
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, wselfall_flag::Bool, prebuilt_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    if prebuilt_flag
        error("Cannot pass prebuilt_flag = true without also passing prebuilt arrays.")
    else
        return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, 
        bnorm_flag, switch_flag, wselfall_flag, true, 
        initialize_prebuilt_arrays(twojmax, length(elements), chem_flag, bzero_flag, bnorm_flag))
    end
end


function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
        chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
        switch_flag::Bool, wselfall_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
        if chem_flag != bnorm_flag
            error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
            bnorm_flag = chem_flag
        end
        return SNAPParams{S, T}(n_atoms, twojmax, elements, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, false)
end
function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T}, chem_flag::Bool,
    bzero_flag::Bool, bnorm_flag::Bool, switch_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAPParams(n_atoms, twojmax, elements, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
    false, false, false)
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAPParams(n_atoms, twojmax, elements, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
    false, false, false)
end
function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool) where {S<:Integer,T<:AbstractFloat}

    # setting bzero_flag == chem_flag
    bzero_flag = chem_flag
    return SNAPParams(n_atoms, twojmax, elements, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bzero_flag,
    false, false, false)
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T}, chem_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, false, bnorm_flag,
    false, false, false)
end

function SNAPParams(n_atoms::S, twojmax::S, elements::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},) where {S<:Integer,T<:AbstractFloat}
    return SNAPParams(n_atoms, twojmax, elements, 
    rcutfac, rmin0, rfac0, radii, weight, false, false, false,
    false, false, false)
end
