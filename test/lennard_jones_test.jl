import Potentials

N = 3
r = Vector{Potentials.Position}(undef, N)
for j = 1:N
    θ = 2*π*rand()
    ϕ = π*rand()
    x = cos(θ)*sin(ϕ)
    y = sin(θ)*sin(ϕ)
    z = cos(ϕ)
    r[j] = Potentials.Position(x, y, z)
end

ϵ = 1.0
σ = 0.5
lj = Potentials.LennardJones(ϵ, σ)
println("Configuration")
show(stdout, "text/plain", r)
println(" ")
println("Lennard Jones Energy ", Potentials.potential_energy(r, lj))
println("Lennard Jones Forces")
show(stdout, "text/plain", Potentials.force(r, lj))
