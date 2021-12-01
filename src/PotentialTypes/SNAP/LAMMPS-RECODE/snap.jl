include("utils.jl")
include("wigner_utils.jl")

function compute_sna(snap::SNAPParams)

   runtime_arrays = initialize_runtime_arrays(snap)

    # These must be done with each new configuration
    for i in 1:snap.n_atoms
        compute_ui(i, 1, snap, runtime_arrays)
    end
    compute_zi(snap, runtime_arrays)
    return runtime_arrays
end 

@time snap = SNAPParams(2, 4, 1, 2.25, 0.1, 0.985)
snap.prebuilt_arrays.rij[1, :] = [1.0, 0.0, 0.0]
snap.prebuilt_arrays.rij[2, :] = [0.0, -1.0, 0.0]
snap.prebuilt_arrays.wj[1]     =  1.0
snap.prebuilt_arrays.wj[2]     =  1.0

snap.prebuilt_arrays.rcutij[1]    = 2.25
snap.prebuilt_arrays.rcutij[2]    = 2.25
# snap.prebuilt_arrays.element[1]   = 1
# snap.prebuilt_arrays.element[2]   = 1
@time runtime_arrays = compute_sna(snap)