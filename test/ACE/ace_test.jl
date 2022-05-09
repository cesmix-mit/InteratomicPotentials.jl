using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials
using InteratomicBasisPotentials
using ACE1

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

ace_params = ACEParams([:Ar], 2, 8, 1.0, 1.0, 0.4, 2.0)
@test isa(ace_params, ACEParams)
e = evaluate_basis(system, ace_params)
@test isa(e, AbstractVector)
f = evaluate_basis_d(system, ace_params)
@test isa(e, AbstractVector)

v = evaluate_basis_v(system, ace_params)
@test isa(e, AbstractVector)

B, dB, W = evaluate_full(system, ace_params)

c = ones(length(ace_params))
ace = ACE(c, ace_params)
@test isa(ace, BasisPotential)

basis = get_rpi(ace_params)
IP = JuLIP.MLIPs.combine(basis, c)
ACE1.Export.export_ace("test.ace", IP)



