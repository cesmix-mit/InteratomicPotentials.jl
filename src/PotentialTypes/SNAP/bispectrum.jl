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

function get_bispectrum(c::Configuration, snap::SNAP; dim = 3)
    J = snap.twojmax
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer multiple of the number of atom types!")
    end
    
    radii_string = "" 
    weight_string = ""
    for j = 1:c.num_atom_types
        radii_string *= string(c.radii[j]) 
        radii_string *= " "
        weight_string *= string(c.weights[j])
        weight_string *= " "
    end
    
    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, num_coeffs, c.num_atoms)
        for i = 1
            # Initialize
            command(lmp, "log none")
            command(lmp, "units " * c.units)
            command(lmp, "boundary $(c.boundaries[1]) $(c.boundaries[2]) $(c.boundaries[3])")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")

            # Setup box
            command(lmp, "region mybox block $(c.x_bounds[1]) $(c.x_bounds[2]) $(c.y_bounds[1]) $(c.y_bounds[2]) $(c.z_bounds[1]) $(c.z_bounds[2])")
            command(lmp, "create_box $(c.num_atom_types) mybox")

            # Create atoms
            for j = 1:c.num_atoms
                atom_id = findall(c.atom_names .== c.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c.Positions[j].x) $(c.Positions[j].y) $(c.Positions[j].z)")
            end

            if c.units == "lj"
                command(lmp, "mass 1 1.0")
            else
                for j = 1:c.num_atom_types
                    command(lmp, "mass $j $(c.Masses[j])")
                end
            end

            # Setup Forcefield
            cutoff = snap.rcutfac * maximum(c.radii)
            max_bounds = max(2*cutoff, min( (c.x_bounds[2] - c.x_bounds[1]), (c.y_bounds[2] - c.y_bounds[1]), (c.z_bounds[2] - c.z_bounds[1]) ) )
            command(lmp, "pair_style zero $max_bounds")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute b all sna/atom $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string
            command(lmp, string_command)
            if dim == 2
                command(lmp, "fix lo all enforce2d")
            end
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

function get_dbispectrum(c::Configuration, snap::SNAP; dim = 3)
    J = snap.twojmax
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer multiple of the number of atom types!")
    end
    
    radii_string = "" 
    weight_string = ""
    for j = 1:c.num_atom_types
        radii_string *= string(c.radii[j]) 
        radii_string *= " "
        weight_string *= string(c.weights[j])
        weight_string *= " "
    end
    
    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, 3*c.num_atom_types*num_coeffs, c.num_atoms)
        for i = 1
            # Initialize
            command(lmp, "log none")
            command(lmp, "units " * c.units)
            command(lmp, "boundary $(c.boundaries[1]) $(c.boundaries[2]) $(c.boundaries[3])")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")

            # Setup box
            command(lmp, "region mybox block $(c.x_bounds[1]) $(c.x_bounds[2]) $(c.y_bounds[1]) $(c.y_bounds[2]) $(c.z_bounds[1]) $(c.z_bounds[2])")
            command(lmp, "create_box $(c.num_atom_types) mybox")

            # Create atoms
            for j = 1:c.num_atoms
                atom_id = findall(c.atom_names .== c.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c.Positions[j].x) $(c.Positions[j].y) $(c.Positions[j].z)")
            end

            if c.units == "lj"
                command(lmp, "mass 1 1.0")
            else
                for j = 1:c.num_atom_types
                    command(lmp, "mass $j $(c.Masses[j])")
                end
            end

            # Setup Forcefield
            cutoff = snap.rcutfac * maximum(c.radii)
            max_bounds = max(2*cutoff, min( (c.x_bounds[2] - c.x_bounds[1]), (c.y_bounds[2] - c.y_bounds[1]), (c.z_bounds[2] - c.z_bounds[1]) ) )
            command(lmp, "pair_style zero $max_bounds")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute db all snad/atom $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
            command(lmp, string_command)
            if dim == 2
                command(lmp, "fix lo all enforce2d")
            end
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "db",  LAMMPS.API.LMP_STYLE_ATOM,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            bispectrum[:, :] = bs
            command(lmp, "clear")
        end
        return bispectrum
    end
    return A

end

function get_vbispectrum(c::Configuration, snap::SNAP; dim = 3)
    J = snap.twojmax
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer multiple of the number of atom types!")
    end
    
    radii_string = "" 
    weight_string = ""
    for j = 1:c.num_atom_types
        radii_string *= string(c.radii[j]) 
        radii_string *= " "
        weight_string *= string(c.weights[j])
        weight_string *= " "
    end
    
    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, 6*c.num_atom_types*num_coeffs, c.num_atoms)
        for i = 1
            # Initialize
            command(lmp, "log none")
            command(lmp, "units " * c.units)
            command(lmp, "boundary $(c.boundaries[1]) $(c.boundaries[2]) $(c.boundaries[3])")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")

            # Setup box
            command(lmp, "region mybox block $(c.x_bounds[1]) $(c.x_bounds[2]) $(c.y_bounds[1]) $(c.y_bounds[2]) $(c.z_bounds[1]) $(c.z_bounds[2])")
            command(lmp, "create_box $(c.num_atom_types) mybox")

            # Create atoms
            for j = 1:c.num_atoms
                atom_id = findall(c.atom_names .== c.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c.Positions[j].x) $(c.Positions[j].y) $(c.Positions[j].z)")
            end

            if c.units == "lj"
                command(lmp, "mass 1 1.0")
            else
                for j = 1:c.num_atom_types
                    command(lmp, "mass $j $(c.Masses[j])")
                end
            end

            # Setup Forcefield
            cutoff = snap.rcutfac * maximum(c.radii)
            max_bounds = max(2*cutoff, min( (c.x_bounds[2] - c.x_bounds[1]), (c.y_bounds[2] - c.y_bounds[1]), (c.z_bounds[2] - c.z_bounds[1]) ) )
            command(lmp, "pair_style zero $max_bounds")
            command(lmp, "pair_coeff * *")
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute vb all snav/atom $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
            command(lmp, string_command)
            if dim == 2
                command(lmp, "fix lo all enforce2d")
            end
            command(lmp, "thermo_style custom pe")
            command(lmp, "run 0")

            ## Extract bispectrum
            bs = extract_compute(lmp, "vb",  LAMMPS.API.LMP_STYLE_ATOM,
                                            LAMMPS.API.LMP_TYPE_ARRAY)
            bispectrum[:, :] = bs
            command(lmp, "clear")
        end
        return bispectrum
    end
    return A

end



function get_snap(c::Configuration, snap::SNAP; dim = 3, reference = true)
    J = snap.twojmax
    if J % 2 == 0
        m = J/2 + 1
        num_coeffs = Int( m * (m+1) * (2*m+1) / 6 )
    elseif J % 2 == 1
        m = (J+1)/2
        num_coeffs = Int( m * (m+1) * (m+2) / 3 )
    else
        AssertionError("twojmax must be an integer multiple of the number of atom types!")
    end
    
    radii_string = "" 
    weight_string = ""
    for j = 1:c.num_atom_types
        radii_string *= string(c.radii[j]) 
        radii_string *= " "
        weight_string *= string(c.weights[j])
        weight_string *= " "
    end
    
    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, c.num_atom_types*num_coeffs + 1, 1+3*c.num_atoms+6)
        for i = 1
            # Initialize
            command(lmp, "log none")
            command(lmp, "units " * c.units)
            command(lmp, "boundary $(c.boundaries[1]) $(c.boundaries[2]) $(c.boundaries[3])")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")

            # Setup box
            command(lmp, "region mybox block $(c.x_bounds[1]) $(c.x_bounds[2]) $(c.y_bounds[1]) $(c.y_bounds[2]) $(c.z_bounds[1]) $(c.z_bounds[2])")
            command(lmp, "create_box $(c.num_atom_types) mybox")

            # Create atoms
            for j = 1:c.num_atoms
                atom_id = findall(c.atom_names .== c.Positions[j].type)[1]
                command(lmp, "create_atoms $atom_id single $(c.Positions[j].x) $(c.Positions[j].y) $(c.Positions[j].z)")
            end

            if c.units == "lj"
                command(lmp, "mass 1 1.0")
            else
                for j = 1:c.num_atom_types
                    command(lmp, "mass $j $(c.Masses[j])")
                end
            end
            # command(lmp, read_data_str)

            # Setup Forcefield
            cutoff = snap.rcutfac * maximum(c.radii)
            max_bounds = max(2*cutoff, min( (c.x_bounds[2] - c.x_bounds[1]), (c.y_bounds[2] - c.y_bounds[1]), (c.z_bounds[2] - c.z_bounds[1]) ) )
            if reference
                command(lmp, "pair_style hybrid/overlay zero $max_bounds zbl 2.0 2.5")
                command(lmp, "pair_coeff * * zero")
                command(lmp, "pair_coeff * * zbl 1 1")
            else
                command(lmp, "pair_style zero $max_bounds")
                command(lmp, "pair_coeff * *")
            end
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute snap all snap $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
            command(lmp, string_command)
            if dim == 2
                command(lmp, "fix lo all enforce2d")
            end
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

    return copy(transpose(A))
end

function get_snap(r::Vector{Configuration}, p; dim = 3, reference = true)
    n = length(r)
    l = length(p.β)
    Aenergy = Array{Float64}(undef, 0, l)
    Aforce = Array{Float64}(undef, 0, l)
    Astress = Array{Float64}(undef, 0, l)
    for j = 1:n
       A = get_snap(r[j], p; dim = dim, reference = reference)
       Aenergy = vcat(Aenergy, reshape(A[1, :], 1, l))
       Aforce = vcat(Aforce, A[2:end-6, :])
       Astress = vcat(Astress, A[end-5:end, :])
    end
    A = vcat(Aenergy, Aforce, Astress)
    return A
end