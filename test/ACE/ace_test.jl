using AtomsBase
using StaticArrays
using Unitful
using UnitfulAtomic
using InteratomicPotentials

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

ace = ACE( species = [:Ar],
           body_order = 2,
           polynomial_degree = 8,
           rcutoff = 2.0)
@test isa(ace, BasisSystem)
e_descrs_l = compute_local_descriptors(system, ace)
@test all(isa.(e_descrs_l, Vector{<:Real}))
@test length(e_descrs_l) == 3
@test length(e_descrs_l[1]) == 8
e_descrs_g = sum(e_descrs_l)
@test isa(e_descrs_g, AbstractVector)

f_descrs = compute_force_descriptors(system, ace)
@test typeof(f_descrs) <: Vector{<:Vector{<:Vector{<:Real}}}
@test length(f_descrs) == 3 # num of atoms
@test length(f_descrs[1]) == 3 # x,y,z
@test length(f_descrs[1][1]) == 8 # number of descriptors

v = compute_virial_descriptors(system, ace)
@test isa(v, AbstractArray)

#print(ace)
lbp = LBasisPotential(ace)
@test isa(lbp, BasisPotential)

basis = get_rpi(ace);

@testset "Reference ACE calculations" begin
    include("reference_ace_test.jl")
end
