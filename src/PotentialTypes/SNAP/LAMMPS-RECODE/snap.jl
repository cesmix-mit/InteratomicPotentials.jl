include("types/types.jl")
include("utilities/utl.jl")

function compute_sna(A::AbstractSystem, snap::SNAPParams)
    # B = Vector{AbstractFloat}[]
    num_coeff = get_num_coeffs(snap.twojmax, length(snap.elements))
    B = [zeros(num_coeff) for i = 1:length(A.particles)]
    dB = [zeros(num_coeff, 3) for i = 1:length(A.particles)]
    W = [zeros(num_coeff, 6) for i = 1:length(A.particles)]
    for  (i, ai) in enumerate(A.particles)
        runtime_arrays = initialize_runtime_arrays(i, ai, A, snap)
        # These must be done with each new configuration
        compute_ui(snap, runtime_arrays)
        compute_dui(snap, runtime_arrays)
        println("ui $(runtime_arrays.ulisttot_r)")
        compute_zi(snap, runtime_arrays)
        println("zi $(runtime_arrays.zlist_r)")
        compute_bi(snap, runtime_arrays)

        for ej in runtime_arrays.element
            compute_dbidrj(ej, snap, runtime_arrays)
        end
        # Energies
        B[i] = runtime_arrays.blist

        # Compute Forces and Stresses
        dB[i] += runtime_arrays.dblist 

        xi, yi, zi = ustrip.(ai.position)
        W[i] += reshape([runtime_arrays.dblist[:, 1]*xi; 
                     runtime_arrays.dblist[:, 2]*yi;
                     runtime_arrays.dblist[:, 3]*zi;
                     runtime_arrays.dblist[:, 2]*zi;
                     runtime_arrays.dblist[:, 1]*zi;
                     runtime_arrays.dblist[:, 1]*yi], num_coeff, 6)

        # These need to be adjusted for multielement systems
        # I think this simple set up applies the wrong set of dblist
        for ind in runtime_arrays.indij
            dB[ind[2]] -= runtime_arrays.dblist

            xj, yj, zj = ustrip.(A.particles[ind[2]].position) 
            W[ind[2]] -= reshape([runtime_arrays.dblist[:, 1]*xj; 
                     runtime_arrays.dblist[:, 2]*yj;
                     runtime_arrays.dblist[:, 3]*zj;
                     runtime_arrays.dblist[:, 2]*zj;
                     runtime_arrays.dblist[:, 1]*zj;
                     runtime_arrays.dblist[:, 1]*yj], num_coeff, 6)
        end
    end
    return B, dB, W
end 


