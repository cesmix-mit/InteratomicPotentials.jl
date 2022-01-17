using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

B = zeros(length(0.1:0.1:2.5), 14)
dB = zeros(length(0.1:0.1:2.5), 14)
for (i,d) in enumerate(0.1:0.1:2.5)
    position1 = @SVector [0.0, 0.0, 0.0] 
    position2 = @SVector [0.0, d, 0.0]

    element  = :Ar
    atoms = [AtomsBase.Atom(element, position1 * 1u"Å"), AtomsBase.Atom(element, position2 * 1u"Å")]

    box = [[-4.0, 4.0], [-4.0, 4.0]]
    system   = FlexibleSystem(box * 1u"Å", [Periodic(), Periodic()], atoms)

    snap = SNAPParams(5, 4, [:Ar], [4.0], 0.01, 0.989, [1.0])
    Btemp, dBtemp, W = compute_sna(system, snap)
    # print(Btemp)
    B[i, :] = sum(Btemp)
    dB[i, :] = dBtemp[1][:, 2]
end
show(stdout, "text/plain", B)
println(" ")

show(stdout, "text/plain", dB)

print(dB)

# show(stdout, "text/plain", dB)
# println(" ")

# show(stdout, "text/plain", W)
# println(" ")
