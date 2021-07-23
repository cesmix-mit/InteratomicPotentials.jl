################################################################################
#
#    This file test some GaN implementations and the ability to build 
#       interatomic potentials through this interface.
#    This file also demonstrates how to use these implementations in order to 
#       solve for the biharmonic coefficients for snap.
#
################################################################################
println("Beginning test of GaN and SNAP implementation.")

import Potentials
using LinearAlgebra: svdvals, cond, pinv
using GalacticOptim, Optim

file_path = "DATA"
c = Potentials.Configuration(file_path; atom_names = [:Ga, :N], 
                    rcutoff = 0.5, neighbor_weight = [1.0, 0.5])


# Set up true potential
ϵ_Ga_Ga = 0.643
σ_Ga_Ga = 2.390
ϵ_N_N   = 1.474
σ_N_N   = 1.981
A_Ga_N  = 608.54
ρ_Ga_N  = 0.435
q_Ga    = 3.0
q_N     = -3.0
ϵ0      = 55.26349406

lj_Ga_Ga = Potentials.LennardJones(ϵ_Ga_Ga,σ_Ga_Ga)
lj_N_N   = Potentials.LennardJones(ϵ_N_N, σ_N_N)
bm_Ga_N  = Potentials.BornMayer(A_Ga_N,ρ_Ga_N)
c_Ga_N   = Potentials.Coulomb(q_Ga, q_N, ϵ0)
gan = Potentials.GaN(lj_Ga_Ga, lj_N_N, bm_Ga_N, c_Ga_N)
println(gan)
pe = Potentials.potential_energy(c, gan)
println("GaN Energy ", pe)
fo = Potentials.force(c, gan)
println("GaN Forces ", fo[1:3])
v = Potentials.virial(c, gan)
println("GaN virial ", v)
v_tensor = Potentials.virial_stress(c, gan)
println("GaN Virial Tensor ", v_tensor)
println(" ")
# Form right hand side
b = [pe; vcat(fo...); vec(v_tensor)]
# b = [pe; vec(v_tensor)]
println("Right hand side shape: ", size(b))

# Form SNAP A matrix
rcutoff = 3.5
twojmax = 3
snap = Potentials.SNAP(rcutoff, twojmax, c.num_atom_types)
A = Potentials.get_snap(c, snap)
# A = [A[1, :]'; A[end-5:end, :]]
println("A matrix shape: ", size(A))
println(" ")
println("β vector shape: ", size(snap.β))

# Look at A diagnostics
println("svdvals of A = ", svdvals(A))
println("Condition # = ", cond(A))

# Fit
snap.β = A \ b 
println("Fitted β = ")
show(stdout, "text/plain", snap.β)
println(" ")
# Test 
pe_snap = Potentials.potential_energy(c, snap)
println("SNAP potential_energy = ", pe_snap)
println("Relative Error = ", abs(pe - pe_snap)/abs(pe))
v_snap = Potentials.virial(c, snap)
println("SNAP Virial = ", v_snap)
println("Virial Error = ", abs(v - v_snap)/abs(v))

# Regularize 
println("Regularizing")
AA = A' * A 
println(" ")
snap.β = pinv(AA, 1e-6) * (A'*b)
println("Regularized Fitted β = ")
show(stdout, "text/plain", snap.β)
println(" ")
pe_snap = Potentials.potential_energy(c, snap)
println("Regularized SNAP potential_energy = ", pe_snap)
println("Regularized Relative Error = ", abs(pe - pe_snap)/abs(pe))
v_snap = Potentials.virial(c, snap)
println("Regularized SNAP Virial = ", v_snap)
println("Regularized Virial Error = ", abs(v - v_snap)/abs(v))

println("End of test.")