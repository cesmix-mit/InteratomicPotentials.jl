##################### Configurations ##########################################
mutable struct Configuration
    file_path                   :: String
    units                       :: String 
    boundaries                  :: Vector{String}
    num_atoms                   :: Int
    num_bonds                   :: Int
    num_angles                  :: Int 
    num_dihedrals               :: Int
    num_impropers               :: Int 
    num_atom_types              :: Int
    num_bond_types              :: Int 
    num_angle_types             :: Int
    num_dihedral_types          :: Int
    num_improper_types          :: Int 
    atom_names                  :: Vector{Symbol}
    radii                       :: Vector{Float64}
    weights                     :: Vector{Float64}
    x_bounds                    :: Vector{Float64}
    y_bounds                    :: Vector{Float64}
    z_bounds                    :: Vector{Float64}
    Masses                      :: Vector{Float64}
    NonBond_coeffs              :: Vector{Float64}
    Bond_coeffs                 :: Vector{Float64}
    Angle_coeffs                :: Vector{Float64}
    Dihedral_coeffs             :: Vector{Float64}
    Improper_coeffs             :: Vector{Float64}
    BondBond_coeffs             :: Vector{Float64}  
    BondAngle_coeffs            :: Vector{Float64}  
    MiddleBondTorsion_coeffs    :: Vector{Float64}
    EndBondTorsion_coeffs       :: Vector{Float64}
    AngleTorsion_coeffs         :: Vector{Float64}
    AngleAngleTorsion_coeffs    :: Vector{Float64}
    BondBond13_coeffs           :: Vector{Float64}
    AngleAngle_coeffs           :: Vector{Float64}
    Positions                   :: Vector{Position}
    Velocities                  :: Vector{Position}
    Bonds                       :: Vector{Float64}
    Angles                      :: Vector{Float64}
    Dihedrals                   :: Vector{Float64}
    Impropers                   :: Vector{Float64}
end

function Configuration()
    file_path           = "" #                   :: String
    units               = "" #                   :: String 
    boundaries          = ["", "", ""] #         :: Vector{String}
    num_atoms           = 0  #                   :: Int
    num_bonds           = 0  #                   :: Int
    num_angles          = 0  #                   :: Int 
    num_dihedrals       = 0  #                   :: Int
    num_impropers       = 0  #                   :: Int 
    num_atom_types      = 0  #                   :: Int
    num_bond_types      = 0  #                   :: Int 
    num_angle_types     = 0  #                   :: Int
    num_dihedral_types  = 0  #                   :: Int
    num_improper_types  = 0  #                   :: Int 
    atom_names          = Vector{Symbol}(undef, 0) #           :: Vector{Symbol}
    radii               = Vector{Float64}(undef, 0) #          :: Vector{Float64}
    weights             = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    x_bounds            = Vector{Float64}(undef, 0) #       :: Vector{Float64}
    y_bounds            = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    z_bounds            = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    Masses              = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    NonBond_coeffs      = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    Bond_coeffs         = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    Angle_coeffs        = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    Dihedral_coeffs     = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    Improper_coeffs     = Vector{Float64}(undef, 0) #        :: Vector{Float64}
    BondBond_coeffs     = Vector{Float64}(undef, 0) #        :: Vector{Float64}  
    BondAngle_coeffs    = Vector{Float64}(undef, 0) #        :: Vector{Float64}  
    MiddleBondTorsion_coeffs  = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    EndBondTorsion_coeffs     = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    AngleTorsion_coeffs       = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    AngleAngleTorsion_coeffs  = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    BondBond13_coeffs         = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    AngleAngle_coeffs         = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    Positions                 = Vector{Float64}(undef, 0) #  :: Vector{Position}
    Velocities                = Vector{Float64}(undef, 0) #  :: Vector{Position}
    Bonds                     = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    Angles                    = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    Dihedrals                 = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    Impropers                 = Vector{Float64}(undef, 0) #  :: Vector{Float64}
    return Configuration(file_path,
    units,
    boundaries,
    num_atoms,
    num_bonds,
    num_angles, 
    num_dihedrals,
    num_impropers, 
    num_atom_types,
    num_bond_types,
    num_angle_types,
    num_dihedral_types,
    num_improper_types,
    atom_names,
    radii,
    weights,
    x_bounds,
    y_bounds,
    z_bounds,
    Masses,
    NonBond_coeffs,
    Bond_coeffs,
    Angle_coeffs,
    Dihedral_coeffs,
    Improper_coeffs,
    BondBond_coeffs, 
    BondAngle_coeffs,  
    MiddleBondTorsion_coeffs,
    EndBondTorsion_coeffs,
    AngleTorsion_coeffs,
    AngleAngleTorsion_coeffs,
    BondBond13_coeffs,
    AngleAngle_coeffs,
    Positions,
    Velocities,
    Bonds,
    Angles,
    Dihedrals,
    Impropers)
end


function get_interparticle_distance(c::Configuration)
    return get_interparticle_distance(c.Positions)
end

function get_interparticle_distance(r::Vector{Configuration})
    n = length(r)
    v = Vector{Vector{Float64}}(undef, n)
    for (i, c) in enumerate(r)
        v[i] = get_interparticle_distance(c)
    end
    return v
end

function get_distance(r::Vector{Configuration}; r0 = Position(0.0, 0.0, 0.0) ::Position)
    N = length(r)
    num_atoms = r[1].num_atoms

    pos = Vector{Vector{Float64}}(undef, N)
    for (i, ri) in enumerate(r)
        pos[i] = get_interparticle_distance(ri, r0)
    end
    return pos 
end

function get_speed(r::Vector{Configuration}; r0 = Position(0.0, 0.0, 0.0) ::Position)
    N = length(r)
    num_atoms = r[1].num_atoms

    pos = Vector{Vector{Float64}}(undef, N)
    for (i, ri) in enumerate(r)
        pos[i] = get_interparticle_distance(ri, r0)
    end
    return pos 
end

function extract_positions(r::Configuration)
    return vec(r.Positions)
end

function extract_positions(r::Vector{Configuration})
    pos = Vector{Vector{Float64}}(undef, length(r))
    for (i, c) in enumerate(r)
        pos[i] = extract_positions(c)
    end
    return pos
end

function extract_velocities(r::Configuration)
    return vec(r.Velocities)
end

function extract_velocities(r::Vector{Configuration})
    vel = Vector{Vector{Float64}}(undef, length(r))
    for (i, c) in enumerate(r)
        vel[i] = extract_velocities(c)
    end
    return vel
end