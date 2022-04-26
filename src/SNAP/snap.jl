include("types/types.jl")
include("utilities/utl.jl")

struct SNAP <: BasisPotential
    coefficients :: AbstractVector
    basis_params :: SNAPParams
end

# get_parameters(snap::SNAP) = Parameter{(:β,)}((snap.coefficients,))
# serialize_parameters(p::Vector, snap::SNAP) = Parameter{(:β,)}((p,))
# deserialize_parameters(p::Parameter{(:β,)}, snap::SNAP) = p.β

# get_hyperparameters(snap::SNAP) = Parameter{(:TwoJMax, :Rcutfac, :Rmin0, :Rfac0, :Radii, :Weights, :Chem, :Bzero, :Bnorm, :Switch, :Wself) }((
#     snap.basis_params.twojmax, snap.basis_params.rcutfac, snap.basis_params.rmin0, snap.basis_params.rfac0, snap.basis_params.radii,
#     snap.basis_params.weight, snap.basis_params.chem_flag, snap.basis_params.bzero_flag, snap.basis_params.bnorm_flag,
#     snap.basis_params.switch_flag, snap.basis_params.wselfall_flag))
# set_hyperparameters(p::Parameter{((:TwoJMax, :Rcutfac, :Rmin0, :Rfac0, 
#                                      :Radii, :Weights, :Chem, :Bzero, 
#                                      :Bnorm, :Switch, :Wself))}, 
#                         snap::SNAP) = SNAP(snap.coefficients, SNAPParams(snap.basis_params.n_atoms, 
#                             p.TwoJMax, snap.basis_params.species, p.Rcutfac, p.Rmin0, p.Rfac0,
#                             p.Radii, p.Weights, p.Chem, p.Bzero, p.Bnorm, p.Switch, p.Wself))
# serialize_hyperparameters(p::Vector, snap::SNAP) = Parameter{(:TwoJMax, :Rcutfac, :Rmin0, :Rfac0, :Radii, :Weights, :Chem, :Bzero, :Bnorm, :Switch, :Wself) }((
#     p[1], p[2], p[3], p[4], p[5:5+snap.basis_params.n_atoms],  
#     p[5+snap.basis_params.n_atoms:5+2*snap.basis_params.n_atoms], p[end-4], p[end-3], p[end-2],
#     p[end-1], p[end]))
# deserialize_hyperparameters(p::Parameter{((:TwoJMax, :Rcutfac, :Rmin0, :Rfac0, 
#                             :Radii, :Weights, :Chem, :Bzero,  
#                             :Bnorm, :Switch, :Wself))}, 
#             snap::SNAP) = [p.TwoJMax, p.Rcutfac, p.Rmin0, p.Rfac0, p.Radii..., p.Weights..., p.Chem, p.Bzero, p.Bnorm, p.Switch, p.Wself]


                                                        

function evaluate_basis(A::AbstractSystem, snap::SNAPParams)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = zeros(num_coeff) 

    for  (i, ai) in enumerate(A)
        i_element = findall(x->x==Symbol(atomic_symbol(ai)), snap.species)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        B += runtime_arrays.blist
    end
    return B
end

function evaluate_basis_d(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    dB = [zeros(num_coeff*length(snap.species), 3) for i = 1:number_of_particles]
    W = [zeros(num_coeff*length(snap.species), 6) for i = 1:number_of_particles]

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

            dB[ii][(i_offset+1):(num_coeff+i_offset), :] += runtime_arrays.dblist
            dB[jj][(i_offset+1):(num_coeff+i_offset), :] -= runtime_arrays.dblist

        end
    end
    return dB
end

function evaluate_basis_v(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    W = [zeros(num_coeff*length(snap.species), 6) for i = 1:number_of_particles]

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

            dB[ii][(i_offset+1):(num_coeff+i_offset), :] += runtime_arrays.dblist
            dB[jj][(i_offset+1):(num_coeff+i_offset), :] -= runtime_arrays.dblist

            xi, yi, zi = ustrip.(A.particles[ii].position)
            xj, yj, zj = ustrip.(A.particles[jj].position)

            W[ii][(i_offset+1):(num_coeff+i_offset), :] += reshape([runtime_arrays.dblist[:, 1]*xi; 
                         runtime_arrays.dblist[:, 2]*yi;
                         runtime_arrays.dblist[:, 3]*zi;
                         runtime_arrays.dblist[:, 2]*zi;
                         runtime_arrays.dblist[:, 1]*zi;
                         runtime_arrays.dblist[:, 1]*yi], num_coeff, 6)
            W[jj][(i_offset+1):(num_coeff+i_offset), :] -= reshape([runtime_arrays.dblist[:, 1]*xj; 
                         runtime_arrays.dblist[:, 2]*yj;
                         runtime_arrays.dblist[:, 3]*zj;
                         runtime_arrays.dblist[:, 2]*zj;
                         runtime_arrays.dblist[:, 1]*zj;
                         runtime_arrays.dblist[:, 1]*yj], num_coeff, 6)
        end
    end
    return sum(W)
end

function evaluate_full(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_snap_coeffs(snap.twojmax, length(snap.species), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = [zeros(num_coeff) for i = 1:number_of_particles]
    dB = [zeros(num_coeff*length(snap.species), 3) for i = 1:number_of_particles]
    W = [zeros(num_coeff*length(snap.species), 6) for i = 1:number_of_particles]

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
            compute_duidrj(rij, wj, rcut, jj, snap, runtime_arrays)
            compute_dbidrj(j_element, snap, runtime_arrays)

            dB[ii][(i_offset+1):(num_coeff+i_offset), :] += runtime_arrays.dblist
            dB[jj][(i_offset+1):(num_coeff+i_offset), :] -= runtime_arrays.dblist

            xi, yi, zi = ustrip.(position(A[ii]))
            xj, yj, zj = ustrip.(position(A[jj]))

            W[ii][(i_offset+1):(num_coeff+i_offset), :] += reshape([runtime_arrays.dblist[:, 1]*xi; 
                         runtime_arrays.dblist[:, 2]*yi;
                         runtime_arrays.dblist[:, 3]*zi;
                         runtime_arrays.dblist[:, 2]*zi;
                         runtime_arrays.dblist[:, 1]*zi;
                         runtime_arrays.dblist[:, 1]*yi], num_coeff, 6)
            W[jj][(i_offset+1):(num_coeff+i_offset), :] -= reshape([runtime_arrays.dblist[:, 1]*xj; 
                         runtime_arrays.dblist[:, 2]*yj;
                         runtime_arrays.dblist[:, 3]*zj;
                         runtime_arrays.dblist[:, 2]*zj;
                         runtime_arrays.dblist[:, 1]*zj;
                         runtime_arrays.dblist[:, 1]*yj], num_coeff, 6)
        end
    end
    return B, dB, W
end
