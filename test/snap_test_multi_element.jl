using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position1 = @SVector [0.0, 0.0, 0.0] 
position2 = @SVector [0.5, .50, 0.0]
position3 = @SVector [0.0, -2.0, 0.0]


atoms = [StaticAtom(position1 * 1u"Å", :Ar), StaticAtom(position2 * 1u"Å", :Ar)]#, StaticAtom(position3 * 1u"Å", :Xe)]

box = [[-4.0, 4.0], [-4.0, 4.0]]
system   = FlexibleSystem(box * 1u"Å", [Periodic(), Periodic()], atoms)

snap = SNAPParams(2, 4, [:Ar], [1.5], 0.00, 0.989, [1.0])
B, dB, W = compute_sna(system, snap)


show(stdout, "text/plain", B)
println(" ")

# show(stdout, "text/plain", dB[1][:, 2]')
# println(" ")

# show(stdout, "text/plain", dB[2][:, 2]')
# println(" ")

# show(stdout, "text/plain", dB[3][:, 2]')
# println(" ")
# show(stdout, "text/plain", dB)
# println(" ")

# show(stdout, "text/plain", W)
# println(" ")
