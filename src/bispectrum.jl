################################################################################
#
#    This file contains methods to calculate the bispectrum of a given Configuration  
#       Requires Configuration (see utils.jl)
################################################################################

function get_bispectrum(c::Configuration, snap::SNAP)
    read_data_str = "read_data " * c.file_path
    rcut_factor = snap.rcut_factor 
    twojmax = snap.twojmax
    rcut_string = "" 
    neighbor_weight_string = ""
    for j = 1:c.num_atom_types
        rcut_string *= string(c.rcutoffs[j]) 
        rcut_string *= " "
        neighbor_weight_string *= string(c.neighbor_weights[j])
        neighbor_weight_string *= " "
    end


    command(lmp, "log none")
    command(lmp, "units metal")
    command(lmp, "boundary p p p")
    command(lmp, "atom_style atomic")
    command(lmp, "atom_modify map array")
    command(lmp, read_data_str)
    command(lmp, "pair_style zero $rcut")
    command(lmp, "pair_coeff * *")
    command(lmp, "compute PE all pe")
    command(lmp, "compute S all pressure thermo_temp")
    compute(lmp, "compute b all sna/atom $rcut 0.99363 $twojmax " * rcut_string * neighbor_weight_string * "rmin 0.0 bzeroflag 0 quadraticflag 0 switchflag 1")
    # command(lmp, "compute SNA all sna/atom $rcut 0.99363 $twojmax " * rcut_string * neighbor_weight_string * "rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1")
    # command(lmp, "compute SNAD all snad/atom $rcut 0.99363 $twojmax * rcut_string * neighbor_weight_string * rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1")
    # command(lmp, "compute SNAV all snav/atom $rcut 0.99363 $twojmax * rcut_string * neighbor_weight_string * rmin0 0.0 bzeroflag 0 quadraticflag 0 switchflag 1")
    command(lmp, "thermo_style custom pe")
    command(lmp, "run 0")

    ## Extract bispectrum
    bs = extract_compute(lmp, "SNA",  LAMMPS.API.LMP_STYLE_ATOM,
                                      LAMMPS.API.LMP_TYPE_ARRAY)
    # deriv_bs = extract_compute(lmp, "SNAD", LAMMPS.API.LMP_STYLE_ATOM,
    #                                         LAMMPS.API.LMP_TYPE_ARRAY)
    return bs, deriv_bs

end