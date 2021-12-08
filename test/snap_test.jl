using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

position1 = @SVector [0.0, 0.0, 0.0] 
position2 = @SVector [0.0, 0.10, 0.0]
position3 = @SVector [0.0, 0.20, 0.0]
position4 = @SVector [0.0, 0.30, 0.0]
position5 = @SVector [0.0, 0.40, 0.0]
element  = :Ar
atoms = []
for (ele, pos) in zip( [:Ar, :Ar, :Ar, :Ar, :Ar], [position1, position2, position3, position4, position5] )
    push!(atoms, StaticAtom(pos * 1u"Å", element))
end

box = [[-4.0, 4.0], [-4.0, 4.0]]
system   = FlexibleSystem(box * 1u"Å", [Periodic(), Periodic()], atoms)

snap = SNAPParams(5, 4, [:Ar], [4.0], 0.05, 0.989, [1.0])
B, dB = compute_sna(system, snap)

show(stdout, "text/plain", B)
println(" ")

show(stdout, "text/plain", dB)
println(" ")
