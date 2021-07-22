################################################################################
#
#    This file exhibits a test optimization problem featuring the Lennard Jones 
#       interatomic potential with 38 atoms.
#    The loss function is modeled as 
#       Total_Error = 0.5*Energy_Error^2 + 0.5*Sum ( Force_Error.^2 )
#
#    There are many pieces of this test that are not representative of a proper
#     workflow that this package aims to provide. That is not the point of this test.
#    The point of this test is to make sure the energies, forces, and gradients are 
#     properly implemented.
#
################################################################################
import Potentials
using GalacticOptim, Optim
println("Beginning test of correctness of LJ implementation.")
# Set up Atoms (using a specific configuration)
N = 38
r = Vector{Potentials.Position}(undef, N)
r0 = [0.163946667699999993 0.319538517000000022 1.71226674359999986;
1.18570812759999988 -1.12890812809999996 -0.616849616600000039;
1.47183871269999988 -0.0530362826999999995 0.944308599800000037;
-0.156601128000000006 1.44806844500000009 -0.969234048100000023;
1.43644861639999988 0.381943172100000006 -0.92280851779999995;
-0.659303829199999969 -1.58099572070000005 -0.355825472500000017;
-1.43644861639999988 -0.381943172100000006 0.92280851779999995;
0.659303829199999969 1.58099572070000005 0.355825472600000026;
-1.18570812749999988 1.12890812809999996 0.616849616600000039;
0.156601128099999987 -1.44806844500000009 0.969234048100000023;
-1.47183871259999988 0.0530362826999999995 -0.944308599700000029;
-0.692522445699999989 0.622174532900000021 -1.48126712230000002;
-0.903125893300000038 0.439946115799999982 1.43235882759999988;
-0.943262934499999983 -0.888676767300000003 -1.17530822099999988;
-0.163946667699999993 -0.319538517000000022 -1.71226674359999986;
0.402594671899999978 -1.70081948270000005 -0.0772747840999999935;
-1.71903085060000005 0.307018840199999976 0.107299306900000002;
0.905297373100000025 1.32824468299999987 -0.690683359700000055;
-0.402594671899999978 1.70081948270000005 0.0772747840999999935;
1.71903085060000005 -0.307018840199999976 -0.107299306799999994;
-0.905297373000000016 -1.32824468299999987 0.690683359700000055;
0.943262934599999991 0.888676767300000003 1.17530822109999988;
0.692522445799999997 -0.622174532900000021 1.48126712230000002;
0.903125893400000046 -0.439946115799999982 -1.43235882749999988;
-1.18439743100000006 -0.634064735699999993 -0.125267553799999987;
0.120331581600000001 -1.00573849309999996 -0.891368452299999947;
0.934263338099999952 -0.873132691499999947 0.430486517999999985;
1.18439743110000006 0.634064735699999993 0.125267553799999987;
-0.934263338000000054 0.873132691499999947 -0.430486517899999976;
-0.370465674600000017 -0.501458934099999976 1.19658741639999988;
0.370465674700000025 0.501458934099999976 -1.19658741639999988;
-0.120331581600000001 1.00573849309999996 0.891368452299999947;
0.649793600800000037 -0.185104513500000012 -0.381540884400000002;
-0.405362064899999985 -0.0660416075000000019 -0.658322833699999999;
0.124574169400000004 0.750628853799999995 -0.15200806280000001;
-0.649793600800000037 0.185104513500000012 0.381540884400000002;
0.405362064999999994 0.0660416075000000019 0.658322833699999999;
-0.124574169299999996 -0.750628853799999995 0.15200806280000001]
for j = 1:N
    r[j] = Potentials.Position(r0[j, 1], r0[j, 2], r0[j, 3])
end


# Construct "True" Lennard Jones
ϵ_true = 1.0
σ_true = 1.0
lj = Potentials.LennardJones(ϵ_true, σ_true)

energies = Potentials.potential_energy(r, lj)
forces = Potentials.force(r, lj)
virials = Potentials.virial(r, lj)

# Implement Loss and Gradient of Loss w.r.t. parameters
function loss(ϵ, σ)
    lj_ = Potentials.LennardJones(ϵ, σ)
    loss = 0.5*(Potentials.potential_energy(r, lj_) - energies)^2 
    vec = (Potentials.force(r, lj_) - forces)
    l = length(vec)
    for v in vec
        loss += 0.5*sum( v .^ 2 ) / l
    end

    loss += 0.5*(Potentials.virial(r, lj_) - virials)^2
    return loss
end


function dloss(ϵ, σ)
    lj_ = Potentials.LennardJones(ϵ, σ)
    grad_potential = Potentials.grad_potential_energy(r, lj_)
    grad_forces = Potentials.grad_force(r, lj_)
    grad_virial = Potentials.grad_virial(r, lj_)
    dloss_ϵ = (Potentials.potential_energy(r, lj_) - energies) .* grad_potential[:dpdϵ] 
    dloss_σ = (Potentials.potential_energy(r, lj_) - energies) .* grad_potential[:dpdσ]
    
    vec = Potentials.force(r, lj_) - forces

    l = length(vec)
    for (v, dv) in zip(vec, grad_forces)
        dloss_ϵ +=  sum( v .* dv[:dfdϵ] )/l
        dloss_σ +=  sum( v .* dv[:dfdσ] )/l
    end

    dloss_ϵ += (Potentials.virial(r, lj_) - virials) .* grad_virial[:dvdϵ] 
    dloss_σ += (Potentials.virial(r, lj_) - virials) .* grad_virial[:dvdσ]
    return [dloss_ϵ, dloss_σ]
end


# Set up Optim.jl optimization problem
f(x, p) = loss(x[1], x[2]) 

function g!(G, x) 
    a, b = dloss(x[1], x[2])
    G[1] = a
    G[2] = b
end
lower = [1e-3, 1e-3]
upper = [4.0, 4.0]

initial = [1.2, 0.85]
println("Truth: ", (ϵ =  ϵ_true, σ = σ_true))
println("Loss at truth = ", loss(ϵ_true, σ_true))
println(" ")
#############################
ff = OptimizationFunction(f, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(ff, initial, [], lb=lower, ub=upper)
@time sol = solve(prob, SAMIN())
println("Finite Difference Solver: ")
print(sol)
println(" ")
println("AD Optimum: ", (ϵ = sol.u[1], σ = sol.u[2]))
println("Error: ϵ=", abs(sol.u[1] - ϵ_true), "\t σ=", abs(sol.u[2] - σ_true))
println(" ")

#############################
ff = OptimizationFunction(f, grad = g!)
prob = OptimizationProblem(ff, initial, [], lb=lower, ub=upper)
@time sol = solve(prob, SAMIN())
println("Analytic Gradient Solver: ")
print(sol)
println(" ")
println("Analytic Gradient Optimum: ", (ϵ = sol.u[1], σ = sol.u[2]))
println("Error: ϵ=", abs(sol.u[1] - ϵ_true), "\t σ=", abs(sol.u[2] - σ_true))
println("Loss at solution = ", loss(sol.u[1], sol.u[2]))
println("End test.")



