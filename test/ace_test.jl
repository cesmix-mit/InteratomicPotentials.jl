using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials
using InteratomicBasisPotentials
using ACE1
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

v = evaluate_basis_v(system, rpi_params)
@test isa(e, AbstractVector)

B, dB, W = evaluate_full(system, rpi_params)

basis = get_rpi(rpi_params)
c = ones(length(basis))
IP = JuLIP.MLIPs.combine(basis, c)
ACE1.Export.export_ace("test.ace", IP)



