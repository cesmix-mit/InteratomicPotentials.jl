include("types/types.jl")
include("utilities/utl.jl")

function compute_sna(A::AbstractSystem, snap::SNAPParams)
    B = Vector{AbstractFloat}[]
    dB = []
    for  a in A.particles
        runtime_arrays = initialize_runtime_arrays(a, A, snap)

        # These must be done with each new configuration
        compute_ui(snap, runtime_arrays)
        compute_dui(snap, runtime_arrays)
        compute_zi(snap, runtime_arrays)
        compute_bi(snap, runtime_arrays)
        compute_dbidrj(1, snap, runtime_arrays)
        push!(B, runtime_arrays.blist)
        push!(dB, runtime_arrays.dblist)
    end
    return B, dB
end 


