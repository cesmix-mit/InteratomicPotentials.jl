################################################################################
#
#    This file contains utility functions
#           1. Functions for parameter containers (redefined NamedTuples)
#               1.1 Adding parameter containers
#               1.2 Converting parameter containers into vectors.
#           2. Positions
#           3. Configurations
#               - A Configuration is a type that holds information typically passed to LAMMPS.
#                   - A Configuration is really just a julia wrapper for information contained in DATA files.
#                   - Configuration contains Positions (depending on the number of atoms).
#               - Configurations are particularly important for SNAP implementation, since
#                   SNAP potential, forces, virial, and learning makes use of LAMMPS (and hence should 
#                   use configurations, instead of Positions).
#           4. System
#               - A system is a collection of Configurations. (Primarily used for training SNAP potentials)

################################################################################
import Base.+
import Base.-
import Base.*
####################### Parameters #############################################
# Abstract type
import Base.NamedTuple as Parameter

# Operators
⊕(p::Parameter, q::Parameter) = merge(p, q)
≅(p::Tuple, q::Tuple) = length(p) == length(q) && all([size(pp) == size(qq) for (pp, qq) in zip(p, q)])
≅(x::NamedTuple{N,T}, y::NamedTuple{N2,T2}) where {N,T,N2,T2} =
  all( [key1 == key2 for (key1, key2) in zip(keys(x), keys(y))] ) && (values(x) ≅ values(y))

+(p::Tuple, q::Tuple) = (length(p) == length(q)) ? [v1 + v2 for (v1, v2) in zip(values(p), values(q))] : AssertionError("p and q must have the same length and sizes.")
+(p::Parameter, q::Parameter) = ( p ≅ q ) ? Parameter{Tuple(keys(p))}( values(p) + values(q) ) : AssertionError("p and q must have the same keys and sizes.")

-(p::Tuple, q::Tuple) = (length(p) == length(q)) ? [v1 - v2 for (v1, v2) in zip(values(p), values(q))] : AssertionError("p and q must have the same length and sizes.")
-(p::Parameter, q::Parameter) = ( p ≅ q ) ? Parameter{Tuple(keys(p))}( values(p) - values(q) ) : AssertionError("p and q must have the same keys and sizes.")
   
*(a::Real, p::Parameter) = Parameter{Tuple(keys(p))}( a .* values(p) )
        
        


## Convert parameter (named tuple) to vector for use in fitting
function parameter_to_vec( p :: Parameter ) 
    l = length(p)
    v = Vector{Float64}()
    for (k, vals) in zip(keys(p), p)
        if typeof(vals) <: AbstractFloat
            append!(v, vals)
        else
            append!(v, vec(vals))
        end
    end
    return v
end

## Convert vector back to parameter for use after fitting
function vec_to_parameter!(p::Parameter, v::Vector{Float64})  
    i = 1
    j = 1
    t = Vector{Any}(undef, length(p))
    for val in p
        l = length(val)
        if l == 1
            t[j] =  v[i]
            i += 1
        else
            s = size(val)
            t[j] = reshape( v[i:(i+l-1)], s)
            i += l
        end
        j += 1
    end
    p = Parameter{keys(p)}(t)
    return p
end



####################### Positions #############################################
struct Position 
    x :: Float64
    y :: Float64 
    z :: Float64
    type ::Symbol
end
function Position(x, y, z)
    return Position(x, y, z, :nothing)
end
function norm(r::Position)
    return sqrt(r.x^2 + r.y^2 + r.z^2)
end

+(r1::Position, r2::Position) = Position(r1.x+r2.x, r1.y + r2.y, r1.z + r2.z)
-(r1::Position, r2::Position) = Position(r1.x-r2.x, r1.y-r2.y, r1.z-r2.z)

##################### Configurations ##########################################
struct Configuration
    file_path                   :: String
    num_atoms                   :: Int
    num_bonds                   :: Int
    num_angles                  :: Int 
    num_dihdedrals              :: Int
    num_impropers               :: Int 
    num_atom_types              :: Int
    num_bond_types              :: Int 
    num_angle_types             :: Int 
    r_cutoffs                   :: Vector{Float64}
    neighbor_weights            :: Vector{Float64}
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

# Read Configuration information from DATA file at file_path
function Configuration(file_path :: String; atom_names = nothing, rcutoff = 0.5, neighbor_weight = 0.5)
    lines = readlines(file_path)
    strip(lines[1]) == "LAMMPS DATA file" ? nothing : AssertionError("path should point to LAMMPS DATA file")
    
    num_atoms           = parse(Int, split(lines[3], " ")[1])
    num_bonds           = parse(Int, split(lines[4], " ")[1])
    num_angles          = parse(Int, split(lines[5], " ")[1])
    num_dihedrals       = parse(Int, split(lines[6], " ")[1])
    num_impropers       = parse(Int, split(lines[7], " ")[1])
    num_atom_types      = parse(Int, split(lines[9], " ")[1])
    num_bond_types      = parse(Int, split(lines[10], " ")[1])
    num_angle_types     = parse(Int, split(lines[11], " ")[1])
    
    xbounds_vec = split(lines[13], " ")
    ybounds_vec = split(lines[14], " ")
    zbounds_vec = split(lines[15], " ")
    x_bounds            = [parse(Float64, xbounds_vec[1]), parse(Float64, xbounds_vec[2])]
    y_bounds            = [parse(Float64, ybounds_vec[1]), parse(Float64, ybounds_vec[2])]
    z_bounds            = [parse(Float64, zbounds_vec[1]), parse(Float64, zbounds_vec[2])]

    if length(rcutoff) == num_atom_types
        r_cutoffs = rcutoff
    else
        r_cutoffs = rcutoff * ones(num_atom_types)
    end
    if length(neighbor_weight) == num_atom_types
        neighbor_weights = neighbor_weight 
    else
        neighbor_weights = neighbor_weight .* ones(num_atom_types)
    end

    # Read Masses
    Masses = zeros(num_atom_types)
    line_num_masses = findall(lines .== "Masses")[1] + 1
    for i = 1:num_atom_types
        Masses[i] = parse(Float64, split(lines[line_num_masses + i])[2])
    end

    # Read Positions

    # Parse atom names
    if atom_names == nothing
        atom_names = [Symbol(i) for i = 1:num_atom_types]
    end

    Positions = Vector{Position}(undef, num_atoms)
    line_num_atoms = findall(lines .== "Atoms")[1] + 1
    for i = 1:num_atoms
        info = split(lines[line_num_atoms+i])
        Positions[i] = Position(parse(Float64, info[3]), 
                                parse(Float64, info[4]), 
                                parse(Float64, info[5]), 
                                atom_names[parse(Int, info[2])]
                                )
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
    Velocities                  = Vector{Float64}(undef, 0)
    Bonds                       = Vector{Float64}(undef, 0)
    Angles                      = Vector{Float64}(undef, 0)
    Dihedrals                   = Vector{Float64}(undef, 0)
    Impropers                   = Vector{Float64}(undef, 0)

    return Configuration(file_path, num_atoms, num_bonds, num_angles, num_dihedrals,num_impropers, 
                        num_atom_types, num_bond_types, num_angle_types, 
                        r_cutoffs, neighbor_weights, x_bounds, y_bounds, z_bounds, 
                        Masses, NonBond_coeffs, Bond_coeffs, Angle_coeffs,
                        Dihedral_coeffs, Improper_coeffs, BondBond_coeffs, 
                        BondAngle_coeffs, MiddleBondTorsion_coeffs, EndBondTorsion_coeffs,
                        AngleTorsion_coeffs, AngleAngleTorsion_coeffs, BondBond13_coeffs, AngleAngle_coeffs, 
                        Positions, Velocities, Bonds, Angles, Dihedrals, Impropers)

end


