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
using LAMMPS
function get_bispectrum(num_atoms, num_elements, file::String)
    num_coeff = 14
    A = LMP(["-screen", "none"]) do lmp
        bispectrum = Array{Float64}(undef, num_coeff, num_atoms)
        for i = 1
            # Initialize
            command(lmp, "log none")
            command(lmp, "units metal")
            command(lmp, "boundary f f f")
            command(lmp, "atom_style atomic")
            command(lmp, "atom_modify map array")

            # Setup box
            command(lmp, "region mybox block -4 4 -4 4 -4 4")
            command(lmp, "read_data $(file)")

            # command(lmp, "mass 1 1")

            # Setup Forcefield
            command(lmp, "pair_style zero 6.0")
            command(lmp, "pair_coeff * *")
            
            command(lmp, "compute PE all pe")
            command(lmp, "compute S all pressure thermo_temp")
            string_command = "compute b all sna/atom 1.5 0.989 4 1.0 1.0 rmin0 0.00 bnormflag 0 bzeroflag 0 switchflag 0" 
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

# function get_dbispectrum(c::Configuration, snap::SNAP; dim = 3)
#     num_coeffs = get_num_coeffs(snap.twojmax)
    
#     radii_string = "" 
#     weight_string = ""
#     for j = 1:c.num_atom_types
#         radii_string *= string(c.radii[j]) 
#         radii_string *= " "
#         weight_string *= string(c.weights[j])
#         weight_string *= " "
#     end
    
#     A = LMP(["-screen", "none"]) do lmp
#         bispectrum = Array{Float64}(undef, 3*c.num_atom_types*num_coeffs, length(c.Atoms))
#         for i = 1
#             # Initialize
#             command(lmp, "log none")
#             command(lmp, "units " * c.units)
#             command(lmp, "boundary $(c.domain.bound_type[1]) $(c.domain.bound_type[2]) $(c.domain.bound_type[3])")
#             command(lmp, "atom_style atomic")
#             command(lmp, "atom_modify map array")

#             # Setup box
#             command(lmp, "region mybox block $(c.domain.bounds[1][1]) $(c.domain.bounds[1][2]) $(c.domain.bounds[2][1]) $(c.domain.bounds[2][2]) $(c.domain.bounds[3][1]) $(c.domain.bounds[3][2])")
#             command(lmp, "create_box $(c.num_atom_types) mybox")

#             # Create atoms
#             for j = 1:length(c.Atoms)
#                 atom_id = findall(c.atom_names .== c.Atoms[j].Type)[1]
#                 command(lmp, "create_atoms $atom_id single $(c.Atoms[j].Position[1]) $(c.Atoms[j].Position[2]) $(c.Atoms[j].Position[3])")
#             end

#             if c.units == "lj"
#                 command(lmp, "mass 1 1.0")
#             else
#                 for j = 1:c.num_atom_types
#                     command(lmp, "mass $j $(c.masses[j])")
#                 end
#             end

#             # Setup Forcefield
#             cutoff = snap.rcutfac * maximum(c.radii)
#             max_bounds = max(2*cutoff, min( (c.domain.bounds[1][2] - c.domain.bounds[1][1]), (c.domain.bounds[2][2] - c.domain.bounds[2][1]), (c.domain.bounds[3][2] - c.domain.bounds[3][1]) ) )
#             command(lmp, "pair_style zero $max_bounds")
#             command(lmp, "pair_coeff * *")
#             command(lmp, "compute PE all pe")
#             command(lmp, "compute S all pressure thermo_temp")
#             string_command = "compute db all snad/atom $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
#             args_command = "" 
#             fields = fieldnames(typeof(snap.keywords))
#             for f in fields
#                 name = string(f)
#                 value = getfield(snap.keywords, f)
#                 if value == 0
#                     continue
#                 else
#                     value = string(value)
#                     if (name == "chemflag")
#                         value = string(c.num_atom_types) * " "
#                         for j = 0:(c.num_atom_types-1)
#                             value *= string(j) * " "
#                         end
#                     end
#                     args_command *= name * " " * value * " "
#                 end
#             end       
#             string_command *= args_command
#             command(lmp, strip(string_command))
#             if dim == 2
#                 command(lmp, "fix lo all enforce2d")
#             end
#             command(lmp, "thermo_style custom pe")
#             command(lmp, "run 0")

#             ## Extract bispectrum
#             bs = extract_compute(lmp, "db",  LAMMPS.API.LMP_STYLE_ATOM,
#                                             LAMMPS.API.LMP_TYPE_ARRAY)
#             bispectrum[:, :] = bs
#             command(lmp, "clear")
#         end
#         return bispectrum
#     end
#     return A

# end

# function get_vbispectrum(c::Configuration, snap::SNAP; dim = 3)
#     num_coeffs = get_num_coeffs(snap.twojmax)
    
#     radii_string = "" 
#     weight_string = ""
#     for j = 1:c.num_atom_types
#         radii_string *= string(c.radii[j]) 
#         radii_string *= " "
#         weight_string *= string(c.weights[j])
#         weight_string *= " "
#     end
    
#     A = LMP(["-screen", "none"]) do lmp
#         bispectrum = Array{Float64}(undef, 6*c.num_atom_types*num_coeffs, length(c.Atoms))
#         for i = 1
#             # Initialize
#             command(lmp, "log none")
#             command(lmp, "units " * c.units)
#             command(lmp, "boundary $(c.domain.bound_type[1]) $(c.domain.bound_type[2]) $(c.domain.bound_type[3])")
#             command(lmp, "atom_style atomic")
#             command(lmp, "atom_modify map array")

#             # Setup box
#             command(lmp, "region mybox block $(c.domain.bounds[1][1]) $(c.domain.bounds[1][2]) $(c.domain.bounds[2][1]) $(c.domain.bounds[2][2]) $(c.domain.bounds[3][1]) $(c.domain.bounds[3][2])")
#             command(lmp, "create_box $(c.num_atom_types) mybox")

#             # Create atoms
#             for j = 1:length(c.Atoms)
#                 atom_id = findall(c.atom_names .== c.Atoms[j].Type)[1]
#                 command(lmp, "create_atoms $atom_id single $(c.Atoms[j].Position[1]) $(c.Atoms[j].Position[2]) $(c.Atoms[j].Position[3])")
#             end

#             if c.units == "lj"
#                 command(lmp, "mass 1 1.0")
#             else
#                 for j = 1:c.num_atom_types
#                     command(lmp, "mass $j $(c.masses[j])")
#                 end
#             end

#             # Setup Forcefield
#             cutoff = snap.rcutfac * maximum(c.radii)
#             max_bounds = max(2*cutoff, min( (c.domain.bounds[1][2] - c.domain.bounds[1][1]), (c.domain.bounds[2][2] - c.domain.bounds[2][1]), (c.domain.bounds[3][2] - c.domain.bounds[3][1]) ) )
#             command(lmp, "pair_style zero $max_bounds")
#             command(lmp, "pair_coeff * *")
#             command(lmp, "compute PE all pe")
#             command(lmp, "compute S all pressure thermo_temp")
#             string_command = "compute vb all snav/atom $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
#             args_command = "" 
#             fields = fieldnames(typeof(snap.keywords))
#             for f in fields
#                 name = string(f)
#                 value = getfield(snap.keywords, f)
#                 if value == 0
#                     continue
#                 else
#                     value = string(value)
#                     if (name == "chemflag")
#                         value = string(c.num_atom_types) * " "
#                         for j = 0:(c.num_atom_types-1)
#                             value *= string(j) * " "
#                         end
#                     end
#                     args_command *= name * " " * value * " "
#                 end
#             end       
#             string_command *= args_command
#             command(lmp, strip(string_command))
#             if dim == 2
#                 command(lmp, "fix lo all enforce2d")
#             end
#             command(lmp, "thermo_style custom pe")
#             command(lmp, "run 0")

#             ## Extract bispectrum
#             bs = extract_compute(lmp, "vb",  LAMMPS.API.LMP_STYLE_ATOM,
#                                             LAMMPS.API.LMP_TYPE_ARRAY)
#             bispectrum[:, :] = bs
#             command(lmp, "clear")
#         end
#         return bispectrum
#     end
#     return A

# end



# function get_snap(f::String)
#     J = snap.twojmax
#     num_coeffs = length(snap.β)
    
#     radii_string = "" 
#     weight_string = ""
#     for j = 1:c.num_atom_types
#         radii_string *= string(c.radii[j]) 
#         radii_string *= " "
#         weight_string *= string(c.weights[j])
#         weight_string *= " "
#     end
    
#     A = LMP(["-screen", "none"]) do lmp
#         bispectrum = Array{Float64}(undef, num_coeffs, 1+3*length(c.Atoms)+6)
#         for i = 1
#             # Initialize
#             command(lmp, "log none")
#             command(lmp, "units " * c.units)
#             command(lmp, "boundary $(c.domain.bound_type[1]) $(c.domain.bound_type[2]) $(c.domain.bound_type[3])")
#             command(lmp, "atom_style atomic")
#             command(lmp, "atom_modify map array")

#             # Setup box
#             command(lmp, "region mybox block $(c.domain.bounds[1][1]) $(c.domain.bounds[1][2]) $(c.domain.bounds[2][1]) $(c.domain.bounds[2][2]) $(c.domain.bounds[3][1]) $(c.domain.bounds[3][2])")
#             command(lmp, "create_box $(c.num_atom_types) mybox")

#             # Create atoms
#             for j = 1:length(c.Atoms)
#                 atom_id = findall(c.atom_names .== c.Atoms[j].Type)[1]
#                 command(lmp, "create_atoms $atom_id single $(c.Atoms[j].Position[1]) $(c.Atoms[j].Position[2]) $(c.Atoms[j].Position[3])")
#             end

#             if c.units == "lj"
#                 command(lmp, "mass 1 1.0")
#             else
#                 for j = 1:c.num_atom_types
#                     command(lmp, "mass $j $(c.masses[j])")
#                 end
#             end

#             # Setup Forcefield
#             cutoff = snap.rcutfac * maximum(c.radii)
#             max_bounds = max(2*cutoff, min( (c.domain.bounds[1][2] - c.domain.bounds[1][1]), (c.domain.bounds[2][2] - c.domain.bounds[2][1]), (c.domain.bounds[3][2] - c.domain.bounds[3][1]) ) )
#             if reference
#                 command(lmp, "pair_style hybrid/overlay zero $max_bounds zbl 2.0 2.5")
#                 command(lmp, "pair_coeff * * zero")
#                 command(lmp, "pair_coeff * * zbl 1 1")
#             else
#                 command(lmp, "pair_style zero $max_bounds")
#                 command(lmp, "pair_coeff * *")
#             end
#             command(lmp, "compute PE all pe")
#             command(lmp, "compute S all pressure thermo_temp")
#             string_command = "compute snap all snap $(snap.rcutfac) 0.99363 $(snap.twojmax) " * radii_string * weight_string 
#             args_command = "" 
#             fields = fieldnames(typeof(snap.keywords))
#             for f in fields
#                 name = string(f)
#                 value = getfield(snap.keywords, f)
#                 if value == 0
#                     continue
#                 else
#                     value = string(value)
#                     if (name == "chemflag")
#                         value = string(c.num_atom_types) * " "
#                         for j = 0:(c.num_atom_types-1)
#                             value *= string(j) * " "
#                         end
#                     end
#                     args_command *= name * " " * value * " "
#                 end
#             end       
#             string_command *= args_command
#             command(lmp, strip(string_command))
#             if dim == 2
#                 command(lmp, "fix lo all enforce2d")
#             end
#             command(lmp, "thermo_style custom pe")
#             command(lmp, "run 0")

#             ## Extract bispectrum
#             bs = extract_compute(lmp, "snap",  LAMMPS.API.LMP_STYLE_GLOBAL,
#                                             LAMMPS.API.LMP_TYPE_ARRAY)
#             bispectrum[:, :] = bs
#             bispectrum[end, 1] = length(c.Atoms)
#             command(lmp, "clear")
#         end
#         return bispectrum
#     end

#     return copy(transpose(A))
# end

# function get_snap(r::Vector{Configuration}, p; dim = 3, reference = true)
#     n = length(r)
#     l = length(p.β)
#     Aenergy = Array{Float64}(undef, 0, l)
#     Aforce = Array{Float64}(undef, 0, l)
#     Astress = Array{Float64}(undef, 0, l)
#     for j = 1:n
#        A = get_snap(r[j], p; dim = dim, reference = reference)
#        Aenergy = vcat(Aenergy, reshape(A[1, :], 1, l))
#        Aforce = vcat(Aforce, A[2:end-6, :])
#        Astress = vcat(Astress, A[end-5:end, :])
#     end
#     A = vcat(Aenergy, Aforce, Astress)
#     return A
# end


A = get_bispectrum(2, 1, "./lammps_snap/starting_configuration.lj")