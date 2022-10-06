include("types/types.jl")
include("utilities/utl.jl")

length(snap::SNAP) = get_num_snap_coeffs(snap.twojmax, size(snap.species, 1), snap.chem_flag)

function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
    switch_flag::Bool, wselfall_flag, prebuilt_flag::Bool, prebuilt_arrays::PrebuiltArrays{T} ) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAP{S, T}(n_atoms, twojmax, species, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, true, prebuilt_arrays)
    
end

function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
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
        return SNAP{S, T}(n_atoms, twojmax, species, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, 
        bnorm_flag, switch_flag, wselfall_flag, true, 
        initialize_prebuilt_arrays(twojmax, size(species, 1), chem_flag, bzero_flag, bnorm_flag))
    end
end


function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
        chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool,
        switch_flag::Bool, wselfall_flag::Bool ) where {S<:Integer,T<:AbstractFloat}
        if chem_flag != bnorm_flag
            error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
            bnorm_flag = chem_flag
        end
        return SNAP{S, T}(n_atoms, twojmax, species, 
        rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
        switch_flag, wselfall_flag, false)
end
function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T}, chem_flag::Bool,
    bzero_flag::Bool, bnorm_flag::Bool, switch_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAP(n_atoms, twojmax, species, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
    false, false, false)
end

function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool, bnorm_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    if chem_flag != bnorm_flag
        error("chem_flag and bnorm_flag should be the same, setting bnorm_flag == chem_flag")
        bnorm_flag = chem_flag
    end
    return SNAP(n_atoms, twojmax, species, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bnorm_flag,
    false, false, false)
end
function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},
    chem_flag::Bool, bzero_flag::Bool) where {S<:Integer,T<:AbstractFloat}

    # setting bzero_flag == chem_flag
    bzero_flag = chem_flag
    return SNAP(n_atoms, twojmax, species, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, bzero_flag, bzero_flag,
    false, false, false)
end

function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T}, chem_flag::Bool) where {S<:Integer,T<:AbstractFloat}
    return SNAP(n_atoms, twojmax, species, 
    rcutfac, rmin0, rfac0, radii, weight, chem_flag, false, bnorm_flag,
    false, false, false)
end

function SNAP(n_atoms::S, twojmax::S, species::Vector{Symbol}, 
    rcutfac::T, rmin0::T, rfac0::T, radii::Vector{T}, weight ::Vector{T},) where {S<:Integer,T<:AbstractFloat}
    return SNAP(n_atoms, twojmax, species, 
    rcutfac, rmin0, rfac0, radii, weight, false, false, false,
    false, false, false)
end


get_rcutoff(snap::SNAP) = snap.rcutoff
get_species(snap::SNAP) = snap.species

function compute_local_descriptors(A::AbstractSystem, snap::SNAP)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = [zeros(num_coeff) for i = 1:length(A)]

    for  (i, ai) in enumerate(A)
        i_element = findall(x->x==Symbol(atomic_symbol(ai)), snap.species)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        B[i] = runtime_arrays.blist
    end
    return B
end

function compute_force_descriptors(A::AbstractSystem, snap::SNAP)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    dB = [zeros(3, num_coeff*length(snap.species)) for i = 1:number_of_particles]

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==atomic_symbol(A, i), snap.species)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        i_offset = num_coeff * (i_element - 1)
        # Compute Forces and Stresses
        for ind = 1:length(runtime_arrays.indij)
            ij = runtime_arrays.indij[ind]
            ii = ij[1]
            jj = ij[2]
            rij = runtime_arrays.rij[ind]
            wj = runtime_arrays.wj[ind]
            rcut = runtime_arrays.rcutij[ind]
            j_element = findall(x->x==atomic_symbol(A, jj), snap.species)[1]
            ## Need to zero-out dulist and dblist each time

            # Now compute dulist for (i,j)
            compute_duidrj(rij, wj, rcut, jj, snap, runtime_arrays)
            compute_dbidrj(j_element, snap, runtime_arrays)

            dB[ii][:, (i_offset+1):(num_coeff+i_offset)] += runtime_arrays.dblist'
            dB[jj][:, (i_offset+1):(num_coeff+i_offset)] -= runtime_arrays.dblist'

        end
    end
    return dB
end

function compute_virial_descriptors(A::AbstractSystem, snap::SNAP)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    W = [zeros(6, num_coeff*length(snap.species)) for i = 1:number_of_particles]

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==Symbol(ai.element.symbol), snap.species)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        B[i] = runtime_arrays.blist

        i_offset = num_coeff * (i_element - 1)
        # Compute Forces and Stresses
        for ind = 1:length(runtime_arrays.indij)
            ij = runtime_arrays.indij[ind]
            ii = ij[1]
            jj = ij[2]
            rij = runtime_arrays.rij[ind]
            wj = runtime_arrays.wj[ind]
            rcut = runtime_arrays.rcutij[ind]
            j_element = findall(x->x==Symbol(A.particles[jj].element.symbol), snap.species)[1]
            ## Need to zero-out dulist and dblist each time

            # Now compute dulist for (i,j)
            compute_duidrj(rij, wj, rcut, jj, snap, runtime_arrays)
            compute_dbidrj(j_element, snap, runtime_arrays)

            xi, yi, zi = ustrip.(A.particles[ii].position)
            xj, yj, zj = ustrip.(A.particles[jj].position)

            W[ii][:, (i_offset+1):(num_coeff+i_offset)] += reshape([runtime_arrays.dblist[:, 1]*xi; 
                         runtime_arrays.dblist[:, 2]*yi;
                         runtime_arrays.dblist[:, 3]*zi;
                         runtime_arrays.dblist[:, 2]*zi;
                         runtime_arrays.dblist[:, 1]*zi;
                         runtime_arrays.dblist[:, 1]*yi], num_coeff, 6)
            W[jj][:, (i_offset+1):(num_coeff+i_offset)] -= reshape([runtime_arrays.dblist[:, 1]*xj; 
                         runtime_arrays.dblist[:, 2]*yj;
                         runtime_arrays.dblist[:, 3]*zj;
                         runtime_arrays.dblist[:, 2]*zj;
                         runtime_arrays.dblist[:, 1]*zj;
                         runtime_arrays.dblist[:, 1]*yj], num_coeff, 6)
        end
    end
    return sum(W)
end

function compute_all_descriptors(A::AbstractSystem, snap::SNAP)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = [zeros(num_coeff) for i = 1:number_of_particles]
    dB = [zeros(3, num_coeff*length(snap.species)) for i = 1:number_of_particles]
    W = [zeros(6, num_coeff*length(snap.species)) for i = 1:number_of_particles]
    for  (i, ai) in enumerate(A)
        i_element = findall(x->x==Symbol(ai.atomic_symbol), snap.species)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        B[i] = runtime_arrays.blist

        i_offset = num_coeff * (i_element - 1)
        # Compute Forces and Stresses
        for ind = 1:length(runtime_arrays.indij)
            ij = runtime_arrays.indij[ind]
            ii = ij[1]
            jj = ij[2]
            rij = runtime_arrays.rij[ind]
            wj = runtime_arrays.wj[ind]
            rcut = runtime_arrays.rcutij[ind]
            j_element = findall(x->x==Symbol(A[jj].atomic_symbol), snap.species)[1]
            ## Need to zero-out dulist and dblist each time

            # Now compute dulist for (i,j)
            compute_duidrj(rij, wj, rcut, ind+1, snap, runtime_arrays)
            compute_dbidrj(j_element, snap, runtime_arrays)

            dB[ii][:, (i_offset+1):(num_coeff+i_offset)] += runtime_arrays.dblist'
            dB[jj][:, (i_offset+1):(num_coeff+i_offset)] -= runtime_arrays.dblist'

            xi, yi, zi = ustrip.(position(A[ii]))
            xj, yj, zj = ustrip.(position(A[jj]))

            W[ii][:, (i_offset+1):(num_coeff+i_offset)] += reshape([runtime_arrays.dblist[:, 1]*xi; 
                         runtime_arrays.dblist[:, 2]*yi;
                         runtime_arrays.dblist[:, 3]*zi;
                         runtime_arrays.dblist[:, 2]*zi;
                         runtime_arrays.dblist[:, 1]*zi;
                         runtime_arrays.dblist[:, 1]*yi], num_coeff, 6)'
            W[jj][:, (i_offset+1):(num_coeff+i_offset)] -= reshape([runtime_arrays.dblist[:, 1]*xj; 
                         runtime_arrays.dblist[:, 2]*yj;
                         runtime_arrays.dblist[:, 3]*zj;
                         runtime_arrays.dblist[:, 2]*zj;
                         runtime_arrays.dblist[:, 1]*zj;
                         runtime_arrays.dblist[:, 1]*yj], num_coeff, 6)'
        end
    end
    return B, dB, sum(W)
end
