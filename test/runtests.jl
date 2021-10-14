using Test

@time begin
@testset "Potentials.jl" begin
    include("lennard_jones_test.jl")
    # include("lammps_test.jl")
end
end
