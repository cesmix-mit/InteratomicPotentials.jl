##################### Configurations ##########################################

include("atoms.jl")
include("angles.jl")
include("bonds.jl")
include("dihedrals.jl")
include("impropers.jl")
include("domains.jl")

struct Configuration{num_atoms, num_atom_types, num_angle_types, num_bond_types,
        num_improper_types, num_dihedral_types, T}
    Atoms                       :: SVector{num_atoms, Atom}
    Angles                      :: SVector{num_angle_types, Angle}
    Bonds                       :: SVector{num_bond_types, Bond}
    Impropers                   :: SVector{num_improper_types, Improper}
    Dihredrals                  :: SVector{num_dihedral_types, Dihedral}

    atom_names                  :: SVector{num_atom_types, Union{Symbol, Nothing}}
    masses                      :: SVector{num_atom_types, T}
    radii                       :: SVector{num_atom_types, T}
    weights                     :: SVector{num_atom_types, T}
    domain                      :: Domain
    units                       :: String
end


