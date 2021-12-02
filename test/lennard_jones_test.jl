using InteratomicPotentials 
import InteratomicPotentials as Potentials

N = 2
r = Vector{Atom}(undef, N)
mass = 38.0

# r[1] = Atom(mass, [0.0, 1.0, 0.0], zeros(3), :Xe)
r[1] = Atom(mass, [-1.0, 0.0, 0.0], zeros(3), :Xe)
r[2] = Atom(mass, [1.0, 0.0, 0.0], zeros(3), :Xe)

ϵ = 1.0
σ = 1.0
lj = LennardJones(ϵ, σ)
println("Beginning test of LJ implementation")
println("Configuration")
show(stdout, "text/plain", r)
println(" ")

true_pe = 4.0 * 1.0 * (1/(2^12) - 1/(2^6)) 
@test isapprox(true_pe, potential_energy(r, lj))
println("Lennard Jones Energy ", potential_energy(r, lj))


println("Lennard Jones Forces")
true_force = 4 *ϵ * (6/(2^6)-12/(2^12)) / 2.0
@test isapprox(true_force, force(r, lj)[1][1])
@test isapprox(true_force, -force(r, lj)[2][1])
show(stdout, "text/plain", force(r, lj))
println(" ")
println("Lennard Jones Virial ", virial(r, lj))
println(" ")
println("Lennard Jones Virial Stress Tensor")
show(stdout, "text/plain", virial_stress(r, lj))
println(" ")

println("Gradient of Lennard Jones Energy with respect to  ϵ, σ ", grad_potential_energy(r, lj))
println("Gradient of Lennard Jones Forces with respect to ϵ, σ")
show(stdout, "text/plain", grad_force(r, lj))
println(" ")

println("Gradient of Lennard Jones Virial with respect to ϵ, σ", grad_virial(r, lj))
println("End test")

