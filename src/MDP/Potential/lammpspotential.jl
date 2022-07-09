using LAMMPS

function lammpspotential(x, t, a, b, c, pbc, unitstyle, pair_style, pair_coeff)

    lmp = LMP()

    # handle boundary conditions
    bcs = ""
    for i = 1:3
        if pbc[i] == 1
            bcs = bcs * " p"
        else  
            bcs = bcs * " f"
        end
    end

    command(lmp, "units " * unitstyle)
    command(lmp, "boundary " * bcs)
    command(lmp, "atom_style atomic")    
    command(lmp, "atom_modify map array")
    command(lmp, "box tilt large")

    # create simulation box
    xlo = 0.0; ylo = 0.0; zlo = 0.0;
    xhi = a[1]; yhi = b[2]; zhi = c[3];
    xy = b[1]; xz = c[1]; yz = c[2];    
    if (abs(xy) + abs(xz) + abs(yz)) <= 1e-10
        command(lmp, "region simbox block $xlo $xhi $ylo $yhi $zlo $zhi units box")
    else
        command(lmp, "region simbox prism $xlo $xhi $ylo $yhi $zlo $zhi $xy $xz $yz units box")
    end
    command(lmp, "create_box 1 simbox")
    
    # pair potential 
    command(lmp, "pair_style " * pair_style)
    for i = 1:length(pair_coeff)
        command(lmp, "pair_coeff " * pair_coeff[i]) 
    end

    # set mass 
    ntypes = length(unique(t))
    for i = 1:ntypes
        itype = Int32(i)
        command(lmp, "mass $itype 1.0")
    end

    # create atom types and atom positions
    n = size(x,2)
    id = Int32.(Array(1:(n)))    
    v = 0.0*x
    image = zeros(Int32, n)
    bexpand = 0
    LAMMPS.API.lammps_create_atoms(lmp, n, id, Int32.(t), x, v, image, bexpand)

    # run lammps
    command(lmp, "run 0")

    # potential energy
    etot = LAMMPS.API.lammps_get_thermo(lmp, "pe")
    pe = Float64(etot)*n

    # extract forces
    forces = extract_atom(lmp, "f")
    id = extract_atom(lmp, "id")
    pos = extract_atom(lmp, "x")
    type = extract_atom(lmp, "type")

    return pe, forces, id, pos, type

    # pos = extract_atom(lmp, "x")
    # t = extract_atom(lmp, "type")
    # i = extract_atom(lmp, "id")    
    # nlocal = extract_global(lmp, "nlocal")
    # ntypes = extract_global(lmp, "ntypes")
end
    
    
    
    
    
    