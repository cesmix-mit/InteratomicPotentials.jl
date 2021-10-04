################################################################################
#
#    This file test some LAMMPS.jl implementations and the ability to build 
#       interatomic potentials through this interface.
#
################################################################################
println("Beginning test of LAMMPS.jl implementation.")
using InteratomicPotentials
using LinearAlgebra: svdvals
file_path = "DATA/1/DATA"

c = load_lammps(file_path; atom_names = [:Ga, :N], radii = [0.5, 0.5], weights = [1.0, 0.5], boundary_type = ["p", "p", "p"], units = "metal" )
snap = SNAP(3.5, 4, c)
println("SNAP: ", snap)
b = get_bispectrum(c, snap)
db = get_dbispectrum(c, snap)
vb = get_vbispectrum(c, snap)
snap_A = get_snap(c, snap)

println("Bispectrum: ", size(b), " Mean: ", sum(b[:, 1])/length(b[:,1]) )
show(stdout, "text/plain", b[:, 1])
println(" ")

println("dBispectrum: ", size(db), " Mean: ", sum(db[:, 1])/length(db[:,1]) )
show(stdout, "text/plain", db[1:10, 1])
println(" ")

println("vBispectrum: ", size(vb), " Mean: ", sum(vb[:, 1])/length(vb[:,1]) )
show(stdout, "text/plain", vb[1:10, 1])
println(" ")

println("SNAP A matrix: ", size(snap_A))
A_temp = [snap_A[1:3, 1:end]; snap_A[end-7:end, 1:end]]
show(stdout, "text/plain", A_temp)
println(" ")
println("Singular values: ")
show(stdout, "text/plain", svdvals(snap_A))
println(" ")
println("Potential Energy Test ", potential_energy(c, snap))
println("Force Test")
show(stdout, "text/plain", force(c, snap)[1:10])
println(" ")
println("Virials Test")
show(stdout, "text/plain", virial(c, snap))
println(" ")
println("End of test.")