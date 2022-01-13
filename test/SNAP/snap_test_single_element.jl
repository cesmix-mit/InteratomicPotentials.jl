using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic
using Statistics

num_elements = 1
num_atoms = 3
twojmax   = 4
rcutfac      = 1.2
radii      = [1.5]
rcut0   = 0.989
weight    = [1.0]

chem_flag = false
bzero_flag = false
bnorm_flag = false
num_coeffs = get_num_coeffs(twojmax, num_elements, chem_flag)

print_flag = false

file = "./SNAP/lammps_snap/starting_configuration_single_element.lj"

include("./lammps_snap/bispectrum_functions.jl")
position1 = @SVector [0.0, 0.0, 0.0] 
position2 = @SVector [0.5, .40, 0.30]
position3 = @SVector [0.2, 0.25, 0.1]


atoms = [StaticAtom(position1 * 1u"Å", :Ar), StaticAtom(position2 * 1u"Å", :Ar), StaticAtom(position3 * 1u"Å", :Ar)]

box = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]
system   = FlexibleSystem(box * 1u"Å", [DirichletZero(), DirichletZero(), DirichletZero()], atoms)
                                                          
snap = SNAPParams(num_atoms, twojmax, [:Ar], rcutfac, 0.00, rcut0, radii, weight, chem_flag, bzero_flag, bnorm_flag)

B, dB, vB = compute_snap(system, snap)

if print_flag
    println("B")
    show(stdout, "text/plain", B)
    println(" \n ")

    println("A")
    show(stdout, "text/plain", A)
    println(" \n ")
end

@test mean(A[:, 1] - B[1]) < 1e-5

@test mean(A[:, 2] - B[2]) < 1e-5

@test mean(A[:, 3] - B[3]) < 1e-5

if print_flag
    println("dB")
    show(stdout, "text/plain", dB[1])
    println(" \n ")

    println("dA")
    show(stdout, "text/plain", reshape(dA[:, 1], :, 3))
    println(" \n ")
end

@test mean(dB[1] - reshape(dA[:, 1], :, 3)) < 1e-5
@test mean(dB[2] - reshape(dA[:, 2], :, 3)) < 1e-5
@test mean(dB[3] - reshape(dA[:, 3], :, 3)) < 1e-5


if print_flag
    println("vA")
    show(stdout, "text/plain", reshape(vA[:, 2], :, 6))
    println(" \n ")

    println("vB")
    show(stdout, "text/plain", vB[2])
    println("\n")
end

@test mean(vB[1] - reshape(vA[:, 1], :, 6)) < 1e-5
@test mean(vB[2] - reshape(vA[:, 2], :, 6)) < 1e-5
@test mean(vB[3] - reshape(vA[:, 3], :, 6)) < 1e-5

