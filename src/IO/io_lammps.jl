# Input to Configuration

function load_lammps(file_path :: String; atom_names = nothing, radii = [3.5], weights = [1.0], boundary_type = ["p", "p", "p"], units = "lj")
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

    Atoms = Vector{Atom}(undef, dict["atoms"])
    Angles = Vector{Angle}(undef, dict["angles"])
    Bonds = Vector{Bond}(undef, dict["bonds"])
    Impropers = Vector{Improper}(undef, dict["impropers"])
    Dihedrals = Vector{Dihedral}(undef, dict["dihedrals"])

    num_atoms           = dict["atoms"]
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

    domain = Domain( [x_bounds, y_bounds, z_bounds], boundary_type )

    if length(radii) == num_atom_types
        radii = radii
    else
        radii = radii[1] .* ones(num_atom_types)
    end
    if length(weights) == num_atom_types
        weights = weights
    else
        weights = weights .* ones(num_atom_types)
    end

    # Read Masses
    masses = zeros(num_atom_types)
    line_num_masses = findall(occursin.("Masses", lines))[1] + 1
    for i = 1:num_atom_types
        masses[i] = parse(Float64, split(lines[line_num_masses + i])[2])
    end

    # Read Per Atom data

    # Parse atom names
    if atom_names == nothing
        atom_names = [Symbol(i) for i = 1:num_atom_types]
    end

    # Prep for Atom positions
    line_num_atoms = findall(occursin.("Atoms", lines))[1] + 1

    # Prep for Velocities
    ind_velocities = findall(occursin.("Velocities", lines))
    if ~isempty(ind_velocities)
        line_num_velocity = ind_velocity[1] + 1
    end
    
    # Iterate through each atom
    for i = 1:num_atoms
        info = split(lines[line_num_atoms+i])
        position = [parse(Float64, info[3]), parse(Float64, info[4]), parse(Float64, info[5])]
        
        ind_velocities = findall(occursin.("Velocities", lines))
        if ~isempty(ind_velocities)
            info = split(lines[line_num_velocity+i])
            velocity = [parse(Float64, info[3]), parse(Float64, info[4]), parse(Float64, info[5])]
        else
            velocity = zeros(3)
        end
        Atoms[parse(Int, info[1])] = Atom(masses[parse(Int, info[2])], 
                                          position, 
                                          velocity, 
                                          atom_names[parse(Int, info[2])]
                                )
    end

    c = Configuration(Atoms, num_atom_types, 
                    Angles, num_angle_types,
                    Bonds, num_bond_types, 
                    Impropers, num_improper_types, 
                    Dihedrals, num_dihedral_types,
                    atom_names, masses, radii, weights,
                    domain, units)

    return c
end

function load_lammps(file_paths :: Vector{String}; atom_names = nothing, radii = [3.5], weights = [1.0], boundary_type = ["p", "p", "p"])
    n = length(file_paths)

    r = Vector{Configuration}(undef, n)
    for i = 1:n
        r[i] = load_lammps(file_paths[i]; atom_names = nothing, radii = [3.5], weights = [1.0], boundary_type = ["p", "p", "p"])
    end

    return r
end

# SAVE
function save_lammps_data(c::Configuration; file = "./DATA")
    open(file, "w") do f
        write(f, "LAMMPS Description \n")
        write(f, "\n")
        
        num_atoms       = length(c.Atoms)
        num_bonds       = length(c.Bonds)
        num_angles      = length(c.Angles)
        num_dihedrals   = length(c.Dihredrals)
        num_impropers   = length(c.Impropers)
        #Specify holistic information 
        for (d, text) in zip([num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers], 
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
        for (bound, text) in zip(c.domain.bounds, ["xlo xhi", "ylo yhi", "zlo zhi"])
            lo = bound[1]
            hi = bound[2]
            write(f, "$lo $hi $text\n")
        end

        write(f, "\n")

        #Specify Masses
        write(f, "Masses\n")
        write(f, "\n")
        for i = 1:c.num_atom_types
            mass = c.masses[i]
            write(f, "$i $mass\n")
        end


        write(f, "\n")
        #Specify Positions 
        write(f, "Atoms\n")
        write(f, "\n")
        for i = 1:num_atoms
            pos = c.Atoms[i].Position
            atom_id = findall(c.atom_names .== pos.type)[1]
            x, y, z = pos.x, pos.y, pos.z
            write(f, "$i $atom_id $x $y $z\n")
        end

    end

end