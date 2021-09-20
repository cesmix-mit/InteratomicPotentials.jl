##################### Configurations ##########################################

include("atoms.jl")
include("angles.jl")
include("bonds.jl")
include("dihedrals.jl")
include("impropers.jl")
include("domains.jl")

mutable struct Configuration
    Atoms                       :: Vector{Atom}
    num_atom_types              :: Int

    Angles                      :: Vector{Angle}
    num_angle_types             :: Int

    Bonds                       :: Vector{Bond}
    num_bond_types              :: Int

    Impropers                   :: Vector{Improper}
    num_improper_types          :: Int

    Dihredrals                  :: Vector{Dihedral}
    num_dihedral_types          :: Int 

    atom_names                  :: Vector
    masses                      :: Vector
    radii                       :: Vector
    weights                     :: Vector
    domain                      :: Domain
    units                       :: String
end


