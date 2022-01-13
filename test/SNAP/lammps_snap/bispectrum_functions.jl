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
    num_coeff = num_coeffs 
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
            radii_command = ""
            for r in radii
                radii_command *= "$r "
            end
            weight_command = ""
            for w in weight
                weight_command *= "$w "
            end
            string_command = "compute b all sna/atom $rcutfac $rcut0 $twojmax " 
            string_command *= radii_command
            string_command *= weight_command 
            string_command *= "rmin0 0.00 bnormflag $(Int(bnorm_flag)) bzeroflag $(Int(bzero_flag)) switchflag 0" 
            
            if chem_flag 
                chem_command = "chem $(num_elements) "
                for i = 0:(num_elements-1)
                    chem_command *= "$(i) "
                end
                string_command *= chem_command
            end
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

function get_dbispectrum(num_atoms, num_elements, file::String)
    num_coeff = 3*num_coeffs*num_elements
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
            radii_command = ""
            for r in radii
                radii_command *= "$r "
            end
            weight_command = ""
            for w in weight
                weight_command *= "$w "
            end
            string_command = "compute db all snad/atom $rcutfac $rcut0 $twojmax " 
            string_command *= radii_command
            string_command *= weight_command 
            string_command *= "rmin0 0.00 bnormflag $(Int(bnorm_flag)) bzeroflag $(Int(bzero_flag)) switchflag 0" 
            
            if chem_flag 
                chem_command = "chem $(num_elements) "
                for i = 0:(num_elements-1)
                    chem_command *= "$(i) "
                end
                string_command *= chem_command
            end
            command(lmp, string_command)
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

function get_vbispectrum(num_atoms, num_elements, file::String)
    num_coeff = 6*num_coeffs*num_elements
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
            radii_command = ""
            for r in radii
                radii_command *= "$r "
            end
            weight_command = ""
            for w in weight
                weight_command *= "$w "
            end
            string_command = "compute vb all snav/atom $rcutfac $rcut0 $twojmax " 
            string_command *= radii_command
            string_command *= weight_command 
            string_command *= "rmin0 0.00 bnormflag $(Int(bnorm_flag)) bzeroflag $(Int(bzero_flag)) switchflag 0" 
            
            if chem_flag 
                chem_command = "chem $(num_elements) "
                for i = 0:(num_elements-1)
                    chem_command *= "$(i) "
                end
                string_command *= chem_command
            end
            command(lmp, string_command)
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



A = get_bispectrum(num_atoms, num_elements, file)
dA = get_dbispectrum(num_atoms, num_elements, file)
vA = get_vbispectrum(num_atoms, num_elements, file)