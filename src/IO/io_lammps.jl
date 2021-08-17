# Input to Configuration

function load_lammps_DATA(file_path :: String; atom_names = nothing, radii = [3.5], weights =[1.0], boundaries = ["f", "f", "f"], units = "lj")
    lines = readlines(file_path)
    
    dict = Dict([("atoms", 0), ("bonds", 0), ("angles", 0), ("dihedrals", 0), ("impropers", 0), 
            ("atom types", 0), ("bond types", 0), ("angle types", 0), ("dihedral types", 0), ("improper types", 0)])
    
    for key = keys(dict)
        line_no = findall(occursin.(key, lines))
        if isempty(line_no)
            continue
        else
            dict[key] = parse(Int, split(lines[line_no[1]], " ")[1])
        end
    end
    num_atoms           = dict["atoms"]
    num_bonds           = dict["bonds"]
    num_angles          = dict["angles"]
    num_dihedrals       = dict["dihedrals"]
    num_impropers       = dict["impropers"]
    num_atom_types      = dict["atom types"]
    num_bond_types      = dict["bond types"]
    num_angle_types     = dict["angle types"]
    num_dihedral_types  = dict["dihedral types"]
    num_improper_types  = dict["improper types"]

    line_no = findall(occursin.("xlo", lines))[1]
    xbounds_vec = split(lines[line_no], " "); line_no += 1
    ybounds_vec = split(lines[line_no], " ");        line_no += 1
    zbounds_vec = split(lines[line_no], " ");      line_no += 1
    x_bounds            = [parse(Float64, xbounds_vec[1]), parse(Float64, xbounds_vec[2])]
    y_bounds            = [parse(Float64, ybounds_vec[1]), parse(Float64, ybounds_vec[2])]
    z_bounds            = [parse(Float64, zbounds_vec[1]), parse(Float64, zbounds_vec[2])]

    if length(radii) == num_atom_types
        radii = radii
    else
        r_cutoffs = rcutoff .* ones(num_atom_types)
    end
    if length(weights) == num_atom_types
        weights = weights
    else
        weights = weights .* ones(num_atom_types)
    end

    # Read Masses
    Masses = zeros(num_atom_types)
    line_num_masses = findall(occursin.("Masses", lines))[1] + 1
    for i = 1:num_atom_types
        Masses[i] = parse(Float64, split(lines[line_num_masses + i])[2])
    end

    # Read Positions

    # Parse atom names
    if atom_names == nothing
        atom_names = [Symbol(i) for i = 1:num_atom_types]
    end

    Positions = Vector{Position}(undef, num_atoms)
    line_num_atoms = findall(occursin.("Atoms", lines))[1] + 1
    for i = 1:num_atoms
        info = split(lines[line_num_atoms+i])
        Positions[parse(Int, info[1])] = Position(parse(Float64, info[3]), 
                                parse(Float64, info[4]), 
                                parse(Float64, info[5]), 
                                atom_names[parse(Int, info[2])]
                                )
    end

    Velocities = Vector{Position}(undef, num_atoms)
    ind_velocities = findall(occursin.("Velocities", lines))
    if ~isempty(ind_velocities)
        line_num_atoms = findall(occursin.("Velocities", lines))[1] + 1
        for i = 1:num_atoms
            info = split(lines[line_num_atoms+i])
            Velocities[parse(Int, info[1])] = Position(parse(Float64, info[2]), 
                                    parse(Float64, info[3]), 
                                    parse(Float64, info[4]), 
                                    Positions[parse(Int, info[1])].type
                                    )
        end
    end

    # To do: Implement other keys
    NonBond_coeffs              = Vector{Float64}(undef, 0)
    Bond_coeffs                 = Vector{Float64}(undef, 0)
    Angle_coeffs                = Vector{Float64}(undef, 0)
    Dihedral_coeffs             = Vector{Float64}(undef, 0)
    Improper_coeffs             = Vector{Float64}(undef, 0)
    BondBond_coeffs             = Vector{Float64}(undef, 0)
    BondAngle_coeffs            = Vector{Float64}(undef, 0)  
    MiddleBondTorsion_coeffs    = Vector{Float64}(undef, 0)
    EndBondTorsion_coeffs       = Vector{Float64}(undef, 0)
    AngleTorsion_coeffs         = Vector{Float64}(undef, 0)
    AngleAngleTorsion_coeffs    = Vector{Float64}(undef, 0)
    BondBond13_coeffs           = Vector{Float64}(undef, 0)
    AngleAngle_coeffs           = Vector{Float64}(undef, 0)
    # Velocities                  = Vector{Position}(undef, num_atoms)
    Bonds                       = Vector{Float64}(undef, 0)
    Angles                      = Vector{Float64}(undef, 0)
    Dihedrals                   = Vector{Float64}(undef, 0)
    Impropers                   = Vector{Float64}(undef, 0)

    return Configuration(file_path, units, boundaries, num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, 
                        num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types,
                        atom_names, radii, weights, x_bounds, y_bounds, z_bounds, 
                        Masses, NonBond_coeffs, Bond_coeffs, Angle_coeffs,
                        Dihedral_coeffs, Improper_coeffs, BondBond_coeffs, 
                        BondAngle_coeffs, MiddleBondTorsion_coeffs, EndBondTorsion_coeffs,
                        AngleTorsion_coeffs, AngleAngleTorsion_coeffs, BondBond13_coeffs, AngleAngle_coeffs, 
                        Positions, Velocities, Bonds, Angles, Dihedrals, Impropers)

end

# SAVE
function save_lammps_data(c::Configuration; file = "./DATA")
    open(file, "w") do f
        write(f, "LAMMPS Description \n")
        write(f, "\n")
        
        #Specify holistic information 
        for (d, text) in zip([c.num_atoms, c.num_bonds, c.num_angles, c.num_dihedrals, c.num_impropers], 
                             ["atoms", "bonds", "angles", "dihedrals", "impropers"])
            write(f, "$d $text\n")
        end

        write(f, "\n")

        #Specify type information
        for (d, text) in zip([c.num_atom_types, c.num_bond_types, c.num_angle_types, c.num_dihedral_types, c.num_improper_types],
                            ["atom types", "bond types", "angle types", "dihedral types", "improper types"]) 
            if d == 0
                continue
            else
                write(f, "$d $text\n")
            end
        end

        write(f, "\n")

        #Specify bounds 
        for (bound, text) in zip([c.x_bounds, c.y_bounds, c.z_bounds], ["xlo xhi", "ylo yhi", "zlo zhi"])
            lo = bound[1]
            hi = bound[2]
            write(f, "$lo $hi $text\n")
        end

        write(f, "\n")

        #Specify Masses
        write(f, "Masses\n")
        write(f, "\n")
        for i = 1:c.num_atom_types
            mass = c.Masses[i]
            write(f, "$i $mass\n")
        end


        write(f, "\n")
        #Specify Positions 
        write(f, "Atoms\n")
        write(f, "\n")
        for i = 1:c.num_atoms
            pos = c.Positions[i]
            atom_id = findall(c.atom_names .== pos.type)[1]
            x, y, z = pos.x, pos.y, pos.z
            write(f, "$i $atom_id $x $y $z\n")
        end

    end

end