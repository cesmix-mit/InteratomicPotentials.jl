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
using LinearAlgebra
using GalacticOptim, Optim
using Random
using Statistics
using Hyperopt
num_train = 41
num_test = 20
indexes = randperm(num_train + num_test)
train_i = indexes[1:num_train]
test_i = indexes[num_train+1:num_train+num_test]
r_train = Vector{Potentials.Configuration}(undef, num_train)
r_test = Vector{Potentials.Configuration}(undef, num_test)
radii = 1.5
println("Radii $radii")
    for (i, index) = enumerate(train_i)
        file_path = "DATA/$index/DATA"
        r_train[i] = Potentials.load_lammps_DATA(file_path; atom_names = [:Ga, :N], radii = [radii, radii], weights = [1.0, 0.5], boundaries = ["p", "p", "p"], units = "metal" )
    end

    for (i, index) = enumerate(test_i)
        file_path = "DATA/$index/DATA"
        r_test[i] = Potentials.load_lammps_DATA(file_path; atom_names = [:Ga, :N], radii = [radii, radii], weights = [1.0, 0.5], boundaries = ["p", "p", "p"], units = "metal" )
    end

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

pe_train = Potentials.potential_energy(r_train, gan)
f_train = Potentials.force(r_train, gan)
v_train = Potentials.virial_stress(r_train, gan)

pe_test = Potentials.potential_energy(r_test, gan)
f_test = Potentials.force(r_test, gan)
v_test = Potentials.virial(r_test, gan)

b = vcat(pe_train, reduce(vcat, reduce(vcat, f_train)), reduce(vcat, v_train))
# println("Right hand side shape: ", size(b))

# Form SNAP A matrix
function cost(rcutoff, radii)
    for i = 1:length(r_train)
        r_train[i].radii = [radii, radii]
    end
    twojmax = 6
    snap = Potentials.SNAP(rcutoff, twojmax, r_train[1])
    A = Potentials.get_snap(r_train, snap)
    # A = [A[1, :]'; A[end-5:end, :]]
    # println("A matrix shape: ", size(A))
    # println(" ")
    # println("β vector shape: ", size(snap.β))

    # Look at A diagnostics
    # println("Condition # = ", cond(A))

    # Fit
    snap.β = A \ b 
    # println("Fitted β = ")
    # show(stdout, "text/plain", snap.β)
    # println(" ")
    # Test 
    pe_snap = Potentials.potential_energy(r_test, snap)
    f_snap = Potentials.force(r_test, snap)
    ff = vcat( (f_test - f_snap)... ) 
    fff = vcat( f_test... )
    v_snap = Potentials.virial(r_test, snap)

    e_error = @. abs(pe_test - pe_snap) / abs(pe_test)
    e_corr = 0.5+0.5*cor(pe_test, pe_snap)
    f_error = @. norm(ff)  / norm( fff )
    f_corr = 0.5+0.5*cor( vcat(vcat(f_test...)...), vcat(vcat(f_snap...)...) )
    v_error = @. abs( v_snap - v_test ) / abs( v_test )
    v_corr = 0.5+0.5*cor(v_test, v_snap)
    println(mean(e_error), "\t", e_corr, "\t")
    println(mean(f_error), "\t", f_corr, "\t")
    println(mean(v_error), "\t", v_corr, "\t")
    total_error = mean(e_error) / e_corr + mean(f_error) / f_corr + mean(v_error) / v_corr 
    return total_error
end


ho = @hyperopt for i = 100, 
    radii = 0.75:0.01:2.25,
    rcut = (1.0):0.01:2.0
    println(i, "\t", radii, "\t", rcut, "   \t")
    @show cost(rcut, radii)
end

best_params, min_f = ho.minimizer, ho.minimum
println(best_params, "\t", min_f)

cost(best_params[2], best_params[1])

println("End of test.")