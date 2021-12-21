using LAMMPS

function lammpssnapdescriptors(x, t, a, b, c, pbc, unitstyle, pair_style, pair_coeff, snapparam, elemradius, elemweight)

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
    natom = size(x,2)
    id = Int32.(Array(1:(natom)))    
    v = 0.0*x
    image = zeros(Int32, natom)
    bexpand = 0
    LAMMPS.API.lammps_create_atoms(lmp, natom, id, Int32.(t), x, v, image, bexpand)

    command(lmp, "compute PE all pe")
    command(lmp, "compute S all pressure thermo_temp")

    # snap parammeters 
    ntypes = Int32(snapparam[1]);  
    twojmax = Int32(snapparam[2]);  
    rcutfac = snapparam[3];
    rfac0 = snapparam[4];
    rmin0 = snapparam[5];
    bzeroflag = Int32(snapparam[6]);
    switchflag = Int32(snapparam[7]);
    quadraticflag = Int32(snapparam[8]);
    chemflag = Int32(snapparam[9]);
    bnormflag = Int32(snapparam[10]);
    wselfallflag = Int32(snapparam[11]);   
    flags = " bzeroflag $bzeroflag quadraticflag $quadraticflag switchflag $switchflag bnormflag $bnormflag wselfallflag $wselfallflag"

    radius = ""
    weight = ""
    nelements = length(elemradius)
    for i = 1:nelements
        elemr = elemradius[i]
        elemw = elemweight[i]
        radius = radius * " $elemr"
        weight = weight * " $elemw"
    end    
    keys = " rmin0 $rmin0"
    if (nelements > 1) & (chemflag==1)                
        keys = keys * " chem $nelements"
        for i = 1:nelements
            itype = Int32(i)
            keys = keys * " $itype"
        end
    end    

    # compute snap descriptors
    command(lmp, "compute snap all snap $rcutfac $rfac0 $twojmax" * radius * weight * flags * keys)    

    # run lammps
    command(lmp, "thermo_style custom pe")
    command(lmp, "run 0")

    # all : [num_bispectrum + num_reference] x [num_energy + num_forces + num_virials]    
    all = extract_compute(lmp, "snap", LAMMPS.API.LMP_STYLE_GLOBAL,
                                       LAMMPS.API.LMP_TYPE_ARRAY)        

    # extract bispectrum, derivatives, and virials                                       
    bi = all[1:end-1,1]
    bd = all[1:end-1,2:(1+3*natom)]
    bv = all[1:end-1,(2+3*natom):end];

    num_bispectrum = length(bi) 
    bd = reshape(bd, (num_bispectrum, 3, natom))
    bd = permutedims(bd, (2,3,1))
    bv = permutedims(bv, (2,1))

    # extract reference energy, forces, and stresses
    energy = all[end,1]
    forces = reshape(all[end,2:(1+3*natom)], (3, natom))    
    stresses = all[end,(2+3*natom):end];
    
    return bi, bd, bv, energy, forces, stresses
end
        

    
    
    
    