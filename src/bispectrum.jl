################################################################################
#
#    This file contains methods to calculate the bispectrum of a given Configuration  
#       Requires Configuration (see utils.jl). All of the command functions use LAMMPS.jl
#       
#       The following methods are useful but require fitted snap.β
#       get_bispectrum calculates the bispectrum components using LAMMPS
#       get_dbispectrum calculuates the derivative of the bispectrum components (used for forces)
#       get_vbispectrum calculuates the virial contribution of the bispectrum components (used for virial stress tensor)
#
#       get_snap produces the A matrix in the SNAP potential fitting paradigm
#           A β = b 
#           where β is the coefficients of the bispectrum components to produce energies, forces, stresses.
#           For detailed information, see http://dx.doi.org/10.1016/j.jcp.2014.12.018 (Thompson et al. 2014)
#       get_snap is instrumental in the methods that calculate potential_energy, force, virial for the snap potential.
################################################################################

function get_bispectrum(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.r_cutoff_factor 
    twojmax = snap.twojmax
    num_coeffs =Int( ( length(snap.β) - 1 ) / c.num_atom_types )
    rcut_string = "" 
    neighbor_weight_string = ""
    for j = 1:c.num_atom_types
        rcut_string *= string(c.r_cutoffs[j]) 
        rcut_string *= " "
        neighbor_weight_string *= string(c.neighbor_weights[j])
        neighbor_weight_string *= " "
    end

    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, num_coeffs, c.num_atoms)
        for j = 1
            command(lmp, "log none")
            command(lmp, "units metal")
            command(lmp, "boundary p p p")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            command(lmp, read_data_str)
            command(lmp, "pair_style zero $rcut_factor")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute b all sna/atom $rcut_factor 0.99363 $twojmax " * rcut_string * neighbor_weight_string * "rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1"
            command(lmp, string_command)
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "b",  LAMMPS.API.LMP_STYLE_ATOM,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            bispectrum[:, :] = bs
            command(lmp, "clear")
        end
        return bispectrum
    end
    return A

end

function get_dbispectrum(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.r_cutoff_factor 
    twojmax = snap.twojmax
    num_coeffs = ( length(snap.β) - 1 ) * 3 
    rcut_string = "" 
    neighbor_weight_string = ""
    for j = 1:c.num_atom_types
        rcut_string *= string(c.r_cutoffs[j]) 
        rcut_string *= " "
        neighbor_weight_string *= string(c.neighbor_weights[j])
        neighbor_weight_string *= " "
    end

    A = LMP(["-screen", "none"]) do lmp
        dbispectrum = Array{Float64}(undef, num_coeffs, c.num_atoms)
        for j = 1
            command(lmp, "log none")
            command(lmp, "units metal")
            command(lmp, "boundary p p p")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            command(lmp, read_data_str)
            command(lmp, "pair_style zero $rcut_factor")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute db all snad/atom $rcut_factor 0.99363 $twojmax " * rcut_string * neighbor_weight_string * "rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1"
            command(lmp, string_command)
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "db",  LAMMPS.API.LMP_STYLE_ATOM,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            dbispectrum[:, :] = bs
            command(lmp, "clear")
        end
        return dbispectrum
    end
    return A
end

function get_vbispectrum(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.r_cutoff_factor 
    twojmax = snap.twojmax
    num_coeffs = (length(snap.β) - 1) * 6 
    rcut_string = "" 
    neighbor_weight_string = ""
    for j = 1:c.num_atom_types
        rcut_string *= string(c.r_cutoffs[j]) 
        rcut_string *= " "
        neighbor_weight_string *= string(c.neighbor_weights[j])
        neighbor_weight_string *= " "
    end

    A = LMP(["-screen", "none"]) do lmp
        vbispectrum = Array{Float64}(undef, num_coeffs, c.num_atoms)
        for j = 1
            command(lmp, "log none")
            command(lmp, "units metal")
            command(lmp, "boundary p p p")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            command(lmp, read_data_str)
            command(lmp, "pair_style zero $rcut_factor")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute vb all snav/atom $rcut_factor 0.99363 $twojmax " * rcut_string * neighbor_weight_string * "rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1"
            command(lmp, string_command)
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "vb",  LAMMPS.API.LMP_STYLE_ATOM,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            vbispectrum[:, :] = bs
            command(lmp, "clear")
        end
        return vbispectrum
    end
    return A
end


function get_snap(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.r_cutoff_factor 
    J = snap.twojmax/c.num_atom_types
    num_coeffs = Int(floor((J + 1) * (J + 2) * (( J + (1.5)) / 3. )))
    rcut_string = "" 
    neighbor_weight_string = ""
    for j = 1:c.num_atom_types
        rcut_string *= string(c.r_cutoffs[j]) 
        rcut_string *= " "
        neighbor_weight_string *= string(c.neighbor_weights[j])
        neighbor_weight_string *= " "
    end

    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, c.num_atom_types*num_coeffs + 1, 1+3*c.num_atoms+6)
        for j = 1
            command(lmp, "log none")
            command(lmp, "units metal")
            command(lmp, "boundary p p p")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")
            command(lmp, read_data_str)
            command(lmp, "pair_style zero $rcut_factor")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute snap all snap $rcut_factor 0.99363 $(snap.twojmax) " * rcut_string * neighbor_weight_string * "rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1"
            command(lmp, string_command)
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "snap",  LAMMPS.API.LMP_STYLE_GLOBAL,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            bispectrum[:, :] = bs
            bispectrum[end, 1] = c.num_atoms
            command(lmp, "clear")
        end
        return bispectrum
    end
    return A'
end

function get_snap(r::Vector{Configuration}, p)
    n = length(r)
    l = length(p.β)
    A = Array{Float64}(undef, 0, l)
    for j = 1:n
       A = vcat(A, get_snap(r[j], p))  
    end
    return A
end