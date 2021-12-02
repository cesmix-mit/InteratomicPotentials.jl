################################################################################
#
#    This file contains types for a variety of empirical atomic potentials
#
################################################################################

const D = 3
abstract type ArbitraryPotential end
abstract type EmpiricalPotential <:ArbitraryPotential end
abstract type FittedPotential <:ArbitraryPotential end
abstract type MixedPotential <:ArbitraryPotential end


############################### Empirical Potentials ################################################
include("EmpiricalPotentials/lj.jl")
include("EmpiricalPotentials/bm.jl")
include("EmpiricalPotentials/coulomb.jl")
include("EmpiricalPotentials/virials.jl")
include("EmpiricalPotentials/vectorization.jl")

################################ SNAP ##############################################################
include("SNAP/snap.jl")

################################ GaN ###############################################################
include("GaN/gan.jl")


################################################################################
# Update on InteratomicPotentials.jl ###########################################
################################################################################

# Modifications to InteratomicPotentials.jl
# TODO: Discuss parametric types {D}

#abstract type ArbitraryPotential{D} end
#abstract type EmpiricalPotential{D} <: ArbitraryPotential{D} end

#using Unitful

## TODO: Unitful.Energy, Unitful.Length
#mutable struct LennardJones <: EmpiricalPotential
#    ϵ::Unitful.Energy
#    σ::Unitful.Length
#end

export inter_pot_conf

function inter_pot_conf(atomic_confs::Vector)
    atom_names = nothing
    radii = [3.5]
    weights = [1.0]
    boundary_type = ["p", "p", "p"]
    units = "lj"
    v = zeros(D)
    new_atomic_confs = []
    new_atomic_confs = Vector{Configuration}(undef, length(atomic_confs))
    for (idx,c) in enumerate(atomic_confs)

        atoms = Vector{Atom}(undef, 0)
        angles = Vector{Angle}(undef, 0)
        bonds = Vector{Bond}(undef, 0)
        impropers = Vector{Improper}(undef, 0)
        dihedrals = Vector{Dihedral}(undef, 0)

        num_atoms           = length(c)
        num_atom_types      = length(unique(species.(c)))
        num_bond_types      = 0
        num_angle_types     = 0
        num_dihedral_types  = 0
        num_improper_types  = 0

        x_bounds = [0.0, bounding_box(c)[1][1].val]
        y_bounds = [0.0, bounding_box(c)[2][2].val]
        z_bounds = [0.0, bounding_box(c)[3][3].val]
        
        domain = Domain(SVector{D}([x_bounds, y_bounds, z_bounds]), SVector{D}(boundary_type))

        if length(radii) != num_atom_types
            radii = radii[1] .* ones(num_atom_types)
        end
        if length(weights) != num_atom_types
            weights = weights .* ones(num_atom_types)
        end

        masses = unique(atomic_mass.(c))

        if atom_names == nothing
            atom_names = [Symbol(i) for i = 1:num_atom_types]
        end

        atoms = []
        for (s, p) in zip(species(c), position(c))
            p2 = [p[1].val, p[2].val, p[3].val]
            push!(atoms, Atom(s.atomic_mass.val, p2, v, Symbol(1))) #TODO
        end

        new_c = Configuration(atoms, num_atom_types,
                              angles, num_angle_types,
                              bonds, num_bond_types,
                              impropers, num_improper_types,
                              dihedrals, num_dihedral_types,
                              atom_names, masses,
                              radii, weights,
                              domain, units)
        new_atomic_confs[idx] = new_c
    end
    return new_atomic_confs
end


function potential_energy(r::SVector, p::LennardJones)
    d = p.σ / norm(r)
    return 4.0 * p.ϵ * ( d^12 - d^6 )
end

# TODO: New functions needed in "gen_test_data"

function potential_energy(s::AbstractSystem, p::ArbitraryPotential)
    N = size(s)[1]
    return sum([potential_energy(position(getindex(s,i)) - position(getindex(s,j)), p)
                for i in 1:N for j in i+1:N])
end

function forces(s::AbstractSystem, p::ArbitraryPotential)
    N = size(s)[1]
    return [ SVector{D}(zeros(D) * 1.0u"N") for i in 1:N ]
end

function virial(s::AbstractSystem, p::ArbitraryPotential) 
    return 0.0
end

function virial_stress(s::AbstractSystem, p::ArbitraryPotential)
    return SMatrix{D, D}(zeros(D, D)u"N/m^2") # TODO: check this
end

function SNAP(rcutfac::Float64, twojmax::Int, num_atom_types::Int)
    keywords = SNAPkeywords(0, 0, 0, 0, 0)
    num_coeffs = get_num_coeffs(twojmax)
    return SNAP(zeros(num_atom_types * num_coeffs + 1), rcutfac, twojmax, keywords)
end


