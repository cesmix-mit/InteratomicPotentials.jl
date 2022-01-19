function load_xyz_data(file_path :: String; atom_names = nothing, radii = [3.5], weights = [1.0], masses = [1.0], boundary_type = ["p", "p", "p"], units = "lj")
    # frames = read_frames(file_path)
    # n = length(frames)
    # r = Vector{Configuration}(undef, n)

    # if atom_names == nothing
    #     atom_names = unique(frames[1]["arrays"]["species"])
    # end
    # num_atom_types = length(atom_names)

    # if length(radii) != num_atom_types
    #     radii = radii[1] .* ones(num_atom_types)
    # end
    # if length(weights) != num_atom_types
    #     weights = weights .* ones(num_atom_types)
    # end
    # if length(masses) !=  num_atom_types
    #     masses = masses .* ones(num_atom_types)
    # end

    # for i = 1:n
    #     frame = frames[i]
    #     num_atoms = frame["N_atoms"]
    #     Atoms = Vector{Atom}(undef, num_atoms)
    #     for j = 1:num_atoms
    #         arrays = frame["arrays"]
    #         name = arrays["species"][j]
    #         id = findall(name .== atom_names)[1]
    #         pos = arrays["pos"][:, j]
    #         if "vel" in keys(arrays)
    #             vel = arrays[vel][:, j]
    #         else
    #             vel = zeros(3)
    #         end
    #         a = Atom(masses[id], pos, vel, name)
    #         Atoms[j] = a

    #     c = Configuration(Atoms, num_atom_types, 
    #                       Angles, num_angle_types, 
    #                       Bonds, num_bond_types, 
    #                       Impropers, num_improper_types,
    #                       Dihedrals, num_dihedral_types,
    #                       atom_names, masses, radii, weights, 
    #                       domain, units)
    #     r[i] = c
    #     end
    # end

    println("This is not yet implemented")
end

function save_xyz_data(c::Configuration, file_path::String)
    println("This is not yet implemented")
end
