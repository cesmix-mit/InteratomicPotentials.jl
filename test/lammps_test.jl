################################################################################
#
#    This file test some LAMMPS.jl implementations and the ability to build 
#       interatomic potentials through this interface.
#
################################################################################
println("Beginning test of LAMMPS.jl implementation.")
import Potentials
file_path = "DATA"

c = Potentials.Configuration(file_path; atom_names = [:Ga, :N], 
                    rcutoff = 0.5, neighbor_weight = [1.0, 0.5])

snap = Potentials.SNAP(3.5, 6)

b = Potentials.get_bispectrum(c, snap)
db = Potentials.get_dbispectrum(c, snap)
vb = Potentials.get_vbispectrum(c, snap)
println("Bispectrum: ", size(b), " Mean: ", sum(b[:, 1])/length(b[:,1]) )
show(stdout, "text/plain", b[1:10, 1])
println(" ")

println("dBispectrum: ", size(db), " Mean: ", sum(db[:, 1])/length(db[:,1]) )
show(stdout, "text/plain", db[1:10, 1])
println(" ")

println("vBispectrum: ", size(vb), " Mean: ", sum(vb[:, 1])/length(vb[:,1]) )
show(stdout, "text/plain", vb[1:10, 1])
println(" ")

println("End of test.")