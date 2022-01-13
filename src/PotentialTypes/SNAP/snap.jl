include("types/types.jl")
include("utilities/utl.jl")

function potential_energy(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    println(nnlist)
    # Get number of coefficients
    num_coeff = get_num_coeffs(snap.twojmax, length(snap.elements), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = zeros(num_coeff) 

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==Symbol(ai.element.symbol), snap.elements)[1]

        runtime_arrays = initialize_runtime_arrays(i, nnlist, A, snap)
        # These must be done with each new configuration
        compute_ui(i_element, snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)

        B += runtime_arrays.blist
    end
    return B
end 

function force(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    println(nnlist)
    # Get number of coefficients
    num_coeff = get_num_coeffs(snap.twojmax, length(snap.elements), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    dB = [zeros(num_coeff*length(snap.elements), 3) for i = 1:number_of_particles]
    W = [zeros(num_coeff*length(snap.elements), 6) for i = 1:number_of_particles]

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==Symbol(ai.element.symbol), snap.elements)[1]

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
            j_element = findall(x->x==Symbol(A.particles[jj].element.symbol), snap.elements)[1]
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

function virial_stress(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    println(nnlist)
    # Get number of coefficients
    num_coeff = get_num_coeffs(snap.twojmax, length(snap.elements), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    W = [zeros(num_coeff*length(snap.elements), 6) for i = 1:number_of_particles]

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==Symbol(ai.element.symbol), snap.elements)[1]

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
            j_element = findall(x->x==Symbol(A.particles[jj].element.symbol), snap.elements)[1]
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

function virial(A::AbstractSystem, snap::SNAPParams)
    v = virial(A, snap)
    return sum(v, dims = 2)
end

function compute_snap(A::AbstractSystem, snap::SNAPParams)
    number_of_particles = length(A.particles)
    # Produce NeighborList
    nnlist = neighborlist(A, snap)
    # Get number of coefficients
    num_coeff = get_num_coeffs(snap.twojmax, length(snap.elements), snap.chem_flag)

    # Initialize SNAP Bispectrum, dBispectrum, and Stress arrays
    B = [zeros(num_coeff) for i = 1:number_of_particles]
    dB = [zeros(num_coeff*length(snap.elements), 3) for i = 1:number_of_particles]
    W = [zeros(num_coeff*length(snap.elements), 6) for i = 1:number_of_particles]

    for  (i, ai) in enumerate(A.particles)
        i_element = findall(x->x==Symbol(ai.element.symbol), snap.elements)[1]

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
            j_element = findall(x->x==Symbol(A.particles[jj].element.symbol), snap.elements)[1]
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
    return B, dB, W
end 


