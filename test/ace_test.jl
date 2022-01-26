using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials
# # Atom from AtomsBase
# struct Atom{D, L<:Unitful.Length, V<:Unitful.Velocity, M<:Unitful.Mass}
#     position::SVector{D, L}
#     velocity::SVector{D, V}
#     atomic_symbol::Symbol
#     atomic_number::Int
#     atomic_mass::M
#     data::Dict{Symbol, Any}  # Store arbitrary data about the atom.
# end

# # Atoms from JuLIP
# mutable struct Atoms{T <: AbstractFloat} <: AbstractAtoms{T}
#     X::Vector{JVec{T}}              # positions
#     P::Vector{JVec{T}}              # momenta (or velocities?)
#     M::Vector{T}                    # masses
#     Z::Vector{AtomicNumber}         # atomic numbers
#     cell::JMat{T}                   # cell
#     pbc::JVec{Bool}                 # boundary condition
#     calc::Union{Nothing, AbstractCalculator}
#     dofmgr::DofManager
#     data::Dict{Any,JData{T}}
#  end
r0 = 2.0
position0 = @SVector [1.0, 0.0, 0.0]
element = :Ar
atom1 = Atom(element, position0 * u"Å")
atom2 = Atom(element, 0.25 * position0 * u"Å")
atom3 = Atom(element, 0.5 * position0 * u"Å")
atom4 = Atom(element, 0.75 * position0 * u"Å")
atom5 = Atom(element, 1.1 * position0 * u"Å")
atom6 = Atom(element, 2.0 * position0 * u"Å")
# atoms = [atom1, atom2, atom3, atom4, atom5, atom6]
atoms = [atom3, atom4, atom5]
box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
bcs = [Periodic(), Periodic(), Periodic()]
system = FlexibleSystem(atoms, box * u"Å", bcs)

rpi_params = RPIParams([:Ar], 2, 8, 1.0, 1.0, 0.4, 2.0)
@test isa(rpi_params, RPIParams)
e = evaluate_basis(system, rpi_params)
@test isa(e, AbstractVector)
f = evaluate_basis_d(system, rpi_params)
@test isa(e, AbstractVector)

# @test Atoms( ustrip.(AtomsBase.position(system)),
#              ustrip.(AtomsBase.velocity(system)),
#              ustrip.(AtomsBase.atomic_mass(system)),
#              ustrip.(AtomsBase.atomic_number(system)),
#              ustrip.(AtomsBase.bounding_box(system)),
#              AtomsBase.boundary_conditions(system),
#                 )


# positions = ustrip.(AtomsBase.position(system))
# velocities = ustrip.(AtomsBase.velocity(system))
# masses = ustrip.(AtomsBase.atomic_mass(system))
# atomic_number = AtomicNumber.(AtomsBase.atomic_number(system))
# cell = JMat(vcat(box...))
# pbc = JVec([pbc == Periodic() ? true : false for pbc in AtomsBase.boundary_conditions(system)]...)
# julip_atoms = Atoms(X = positions, P = velocities, M = masses, Z = atomic_number, cell = cell, pbc = pbc)

# # Construct RPI basis of Sparse Spherical Harmonics
# rpi_ = rpi_basis(; 
#       species = [:Ar, :Ne],
#       N = 1,                        # correlation order = body-order - 1
#       maxdeg = 8,                  # polynomial degree
#       D = SparsePSHDegree(; wL=1.5, csp=1.0),
#       r0 = r0,                      # estimate for NN distance
#       rin = 0.65*r0, rcut = 5.0,    # domain for radial basis (cf documentation)
#       pin = 0)

# # Get rpi_descriptors
# descritptors = ACE1.energy(rpi_, julip_atoms) # length(rpi_) Vector{Float64}
# grad_descriptors = ACE1.forces(rpi_, julip_atoms) # length(rpi_) Vector{SVector{3, Float64}}


# # Construct RPI basis of Transformed Jacobi
# rcut = 5.5
# rin = 0.0
# maxn = 5
# species = [:Ar, :Ne]
# D = SparsePSHDegree(wL = 1.0)

# trans = ACE1.PolyTransform(1, r0)
# J = ACE1.OrthPolys.transformed_jacobi(10, trans, rcut, rin)
# P1 = ACE1.RPI.PSH1pBasis(J, maxn, D=D, species = species)
# J_basis_ = RPIBasis(P1, 3, D, maxn)

# descritptors = ACE1.energy(J_basis_, julip_atoms) # length(rpi_) Vector{Float64}
# grad_descriptors = ACE1.forces(J_basis_, julip_atoms) # length(rpi_) Vector{SVector{3, Float64}}