using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials

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

ace = ACE( species = [:Ar], body_order = 2, polynomial_degree = 8, 
           wL = 1.0, csp = 1.0, r0 = 0.4, rcutoff = 2.0)
@test isa(ace, BasisSystem)
e = sum(compute_local_descriptors(system, ace))
@test isa(e, AbstractVector)
f = compute_local_descriptors(system, ace)
@test all(isa.(f, (AbstractVector,)))
v = compute_virial_descriptors(system, ace)
@test isa(v, AbstractArray)

print(ace)
lbp = LBasisPotential(ace)
@test isa(lbp, BasisPotential)

basis = get_rpi(ace);


