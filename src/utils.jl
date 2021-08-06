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
#           4. System (not yet implemented)
#               - A system is a collection of Configurations. (Primarily used for training SNAP potentials)

################################################################################
import Base.+
import Base.-
import Base.*
import Base.copy
import Base.vec
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

+(r1::Position, r2::Position) = (r1.type == r2.type ? Position(r1.x+r2.x, r1.y + r2.y, r1.z + r2.z, r1.type) : AssertionError("Atoms must have the same type for this operation!"))
-(r1::Position, r2::Position) = (r1.type == r2.type ? Position(r1.x-r2.x, r1.y - r2.y, r1.z - r2.z, r1.type) : AssertionError("Atoms must have the same type for this operation!"))
*(a::AbstractFloat, r::Position) = Position(a*r.x, a*r.y, a*r.z, r.type)
*(r::Position, a::AbstractFloat) = Position(a*r.x, a*r.y, a*r.z, r.type)

+(a::AbstractFloat, r::Position) = Position(a + r.x, a + r.y, a + r.z, r.type)
+(a::Vector{Float64}, r::Position) = Position(a[1] + r.x, a[2] + r.y, a[3] + r.z, r.type)
+(r::Position, a::Vector{Float64}) = Position(a[1] + r.x, a[2] + r.y, a[3] + r.z, r.type)
-(a::Vector{Float64}, r::Position) = Position(a[1] - r.x, a[2] - r.y, a[3] - r.z, r.type)

vec(r::Position) = [r.x, r.y, r.z]
vec(r::Vector{Position}) = [[ri.x, ri.y, ri.z] for ri in r]

function get_interparticle_distance(ri::Position, rj::Position)
    return norm(ri - rj)
end
function get_interparticle_distance(r::Vector{Position})
    n = length(r)
    distances = zeros(0)    
    for i = 1:n
        ri = r[i]
        for j = (i+1):n 
            d = get_interparticle_distance(ri, r[j])
            push!(distances, d, d)
        end
    end
    return distances
end
##################### Configurations ##########################################
struct Configuration
    file_path                   :: String
    num_atoms                   :: Int
    num_bonds                   :: Int
    num_angles                  :: Int 
    num_dihedrals              :: Int
    num_impropers               :: Int 
    num_atom_types              :: Int
    num_bond_types              :: Int 
    num_angle_types             :: Int
    num_dihedral_types          :: Int
    num_improper_types          :: Int 
    atom_names                  :: Vector{Symbol}
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
function Configuration(file_path :: String; atom_names = nothing, rcutoff = [0.5], neighbor_weight =[1.0])
    lines = readlines(file_path)
    strip(lines[1]) == "LAMMPS DATA file" ? nothing : AssertionError("path should point to LAMMPS DATA file")
    
    num_atoms           = parse(Int, split(lines[3], " ")[1])
    num_bonds           = parse(Int, split(lines[4], " ")[1])
    num_angles          = parse(Int, split(lines[5], " ")[1])
    num_dihedrals       = parse(Int, split(lines[6], " ")[1])
    num_impropers       = parse(Int, split(lines[7], " ")[1])

 
    line_no = 9
    num_types = zeros(Int, 5)
    for (i, (nums, text)) in enumerate(zip([num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers],
                            ["atom", "bond", "angle", "dihedral", "improper"]))
        if nums == 0
            line_length = length(split(lines[line_no], " "))
            if line_length > 1
                line_no +=1
            else
                continue
            end
        
        else text == split(lines[line_no], " ")[2]
            num_types[i] = parse(Int, split(lines[line_no], " ")[1])
            line_no += 1
        end
    end
    num_atom_types = num_types[1]
    num_bond_types = num_types[2]
    num_angle_types = num_types[3]
    num_dihedral_types = num_types[4]
    num_improper_types = num_types[5]

    line_no += 1
    xbounds_vec = split(lines[line_no], " "); line_no += 1
    ybounds_vec = split(lines[line_no], " ");        line_no += 1
    zbounds_vec = split(lines[line_no], " ");      line_no += 1
    x_bounds            = [parse(Float64, xbounds_vec[1]), parse(Float64, xbounds_vec[2])]
    y_bounds            = [parse(Float64, ybounds_vec[1]), parse(Float64, ybounds_vec[2])]
    z_bounds            = [parse(Float64, zbounds_vec[1]), parse(Float64, zbounds_vec[2])]

    if length(rcutoff) == num_atom_types
        r_cutoffs = rcutoff
    else
        r_cutoffs = rcutoff .* ones(num_atom_types)
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
    Velocities                  = Vector{Position}(undef, num_atoms)
    Bonds                       = Vector{Float64}(undef, 0)
    Angles                      = Vector{Float64}(undef, 0)
    Dihedrals                   = Vector{Float64}(undef, 0)
    Impropers                   = Vector{Float64}(undef, 0)

    return Configuration(file_path, num_atoms, num_bonds, num_angles, num_dihedrals, num_impropers, 
                        num_atom_types, num_bond_types, num_angle_types, num_dihedral_types, num_improper_types,
                        atom_names, r_cutoffs, neighbor_weights, x_bounds, y_bounds, z_bounds, 
                        Masses, NonBond_coeffs, Bond_coeffs, Angle_coeffs,
                        Dihedral_coeffs, Improper_coeffs, BondBond_coeffs, 
                        BondAngle_coeffs, MiddleBondTorsion_coeffs, EndBondTorsion_coeffs,
                        AngleTorsion_coeffs, AngleAngleTorsion_coeffs, BondBond13_coeffs, AngleAngle_coeffs, 
                        Positions, Velocities, Bonds, Angles, Dihedrals, Impropers)

end

function save_as_lammps_data(c::Configuration; file = "./DATA")
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

function get_interparticle_distance(c::Configuration)
    return get_interparticle_distance(c.Positions)
end

function get_positions(r::Vector{Configuration})
    N = length(r)
    num_atoms = r[1].num_atoms

    pos = Vector{Vector{Float64}}(undef, N)
    for (i, ri) in enumerate(r)
        pos[i] = Potentials.norm.(ri.Positions)
    end
    return pos 
end

function get_velocities(r::Vector{Configuration})
    N = length(r)
    num_atoms = r[1].num_atoms

    pos = Vector{Vector{Float64}}(undef, N)
    for (i, ri) in enumerate(r)
        pos[i] = Potentials.norm.(ri.Velocities)
    end
    return pos 
end

