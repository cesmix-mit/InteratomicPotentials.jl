########################################################################
##############################  SNAP  ###################################
#########################################################################

mutable struct SNAPkeywords
    n_elements      :: Int
    switchflag      :: Int 
    bzeroflag       :: Int 
    quadraticflag   :: Int
    chemflag        :: Int 
    bnormflag       :: Int
end

mutable struct SNAP <: FittedPotential
    β :: Vector{Float64}
    rcutfac         :: Float64
    twojmax         :: Int
    keywords        :: SNAPkeywords
end

function get_num_coeffs(twojmax :: Int)
    J = twojmax 
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer!")
    end 
    return num_coeffs
end


function SNAP(rcutfac::Float64, twojmax::Int, c::Configuration)
    keywords = SNAPkeywords(0, 0, 0, 0, 0)
    num_coeffs = get_num_coeffs(twojmax)

    return SNAP( zeros(c.num_atom_types * num_coeffs + 1) , rcutfac, twojmax, keywords)
end

function SNAP(rcutfac::Float64, twojmax::Int, c::Configuration, keywords::SNAPkeywords)
    num_coeffs = get_num_coeffs(twojmax)
    if keywords.quadraticflag == 1
        num_coeffs += num_coeffs * (num_coeffs+1) / 2
    end
    if keywords.chemflag == 1
        num_coeffs = num_coeffs * (c.num_atom_types^3)
        keywords.bnormflag = 1
    else
        num_coeffs = num_coeffs * c.num_atom_types
    end

    return SNAP( zeros(num_coeffs + 1) , rcutfac, twojmax, keywords)
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
        for i = 1:c.num_atom_types
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
include("virial.jl")
include("gradients.jl")
