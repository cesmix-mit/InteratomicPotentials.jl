using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

num_elements = 1
num_atoms = 3
twojmax = 4
rcutfac = 1.2
radii = [1.5]
rcut0 = 0.989
weight = [1.0]

chem_flag = false
bzero_flag = false
bnorm_flag = false
snap = SNAP(num_atoms, twojmax, [:Ar], rcutfac, 0.00, rcut0, radii, weight, chem_flag, bzero_flag, bnorm_flag)

num_coeffs = length(snap)

print_flag = false

file = "./SNAP/lammps_snap/starting_configuration_single_element.lj"

include("./lammps_snap/bispectrum_functions.jl")
position1 = @SVector [0.0, 0.0, 0.0]
position2 = @SVector [0.5, 0.40, 0.30]
position3 = @SVector [0.2, 0.25, 0.1]


atoms = [Atom(:Ar, position1 * u"Å"), Atom(:Ar, position2 * u"Å"), Atom(:Ar, position3 * u"Å")]

box = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]
system = FlexibleSystem(atoms, box * u"Å", [DirichletZero(), DirichletZero(), DirichletZero()])

B = compute_local_descriptors(system, snap)

if print_flag
    println("B")
    show(stdout, "text/plain", B)
    println("\n")
    println("length of B", length(B))
    println("length of B[1]", length(B[1]))
end

dB = compute_force_descriptors(system, snap)
if print_flag
    println("dB")
    show(stdout, "text/plain", dB)
    println("\n")
    println("length of dB", length(dB))
    println("length of dB[1]", size(dB[1]))
end

B, dB, vB = compute_all_descriptors(system, snap)

if print_flag
    println("B")
    show(stdout, "text/plain", B)
    println(" \n ")

    println("A")
    show(stdout, "text/plain", A)
    println(" \n ")
end

@test sum(A[:, 1] - B[1]) < 1e-5

@test sum(A[:, 2] - B[2]) < 1e-5

@test sum(A[:, 3] - B[3]) < 1e-5

if print_flag
    println("dB")
    show(stdout, "text/plain", dB[1])
    println(" \n ")

    println("dA")
    show(stdout, "text/plain", reshape(dA[:, 1], :, 3))
    println(" \n ")
end

@test sum(dB[1]' - reshape(dA[:, 1], :, 3)) < 1e-5
@test sum(dB[2]' - reshape(dA[:, 2], :, 3)) < 1e-5
@test sum(dB[3]' - reshape(dA[:, 3], :, 3)) < 1e-5


if print_flag
    println("vA")
    show(stdout, "text/plain", reshape(vA[:, 2], :, 6))
    println(" \n ")

    println("vB")
    show(stdout, "text/plain", vB[2])
    println("\n")
end
vA = sum(vA[:, i] for i = 1:3)
@test sum(vB' - reshape(vA, :, 6)) < 1e-5

