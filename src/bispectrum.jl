################################################################################
#
#    This file contains methods to calculate the bispectrum of a given Configuration  
#       Requires Configuration (see utils.jl)
################################################################################

function get_bispectrum(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.r_cutoff_factor 
    twojmax = snap.twojmax
    num_coeffs = length(snap.β) - 1
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
    num_coeffs = (length(snap.β) - 1) * 3 * c.num_atom_types
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
    num_coeffs = (length(snap.β) - 1) * 6 * c.num_atom_types
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
