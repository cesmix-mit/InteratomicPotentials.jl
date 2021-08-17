########################################################################
##############################  SNAP  ###################################
#########################################################################

mutable struct SNAP <: FittedPotential
    β :: Vector{Float64}
    rcutfac         :: Float64
    twojmax         :: Int
    num_atom_types  :: Int
end

function SNAP(rcutfac::Float64, twojmax::Int, num_atom_types::Int)
    J = twojmax 
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer multiple of the number of atom types!")
    end 
    return SNAP(ones(num_atom_types*num_coeffs+1), rcutfac, twojmax, num_atom_types)
end

function get_trainable_params(snap::SNAP)
    return (β = snap.β, )
end

function set_trainable_params!(snap::SNAP, params::NamedTuple)
    snap.β = params.β
    return snap
end

function get_nontrainable_params(snap::SNAP)
    return (rcut = snap.rcut, twojmax = snap.twojmax)
end

function create_snap_files(c::Configuration, snap::SNAP, file)
    open(file*".snapcoeff", "w") do f 
        write(f, "#UNITS: $(c.units)\n")
        write(f, "#LAMMPS SNAP COEFFICIENTS\n")
        write(f, "$(c.num_atom_types) $(length(snap.β))\n")
        for i = 1:snap.num_atom_types
            write(f, "$(c.atom_names[i]) $(c.radii[i]) $(c.weights[i])\n")
        end
        write(f, "$(snap.β[end])\n")
        for b = snap.β[1:end-1]
            write(f, "$b\n")
        end
    end

    open(file*".snapparam", "w") do f
        write(f, "#UNITS: $(c.units)\n")
        write(f, "#LAMMPS SNAP PARAMETERS\n")
        write(f, "rcutfac $(snap.rcutfac)\n")
        write(f, "twojmax $(snap.twojmax)\n")
    end
end

##########################################################################################

include("bispectrum.jl")
include("energy.jl")
include("force.jl")
include("virial")
