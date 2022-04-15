@testset "Neighbor List Unit Tests" begin
    atoms = [
        Atom(:Ar, (@SVector [0.5, 0.1, 0.3])u"bohr"),
        Atom(:Ar, (@SVector [0.7, 2.5, -2.1])u"bohr"),
        Atom(:Ar, (@SVector [1.1, -3.1, 7.3])u"bohr")
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    boundary_conditions = [Periodic(), DirichletZero(), Periodic()]
    system = FlexibleSystem(atoms, box, boundary_conditions)

    rcutoff = 5.0
    nnlist = neighborlist(system, rcutoff)

    true_j = [
        Int[2, 3],
        Int[],
        Int[]
    ]

    true_R = [
        Float64[norm([-0.2, -2.4, 0.4]), norm([0.4, 3.2, 0.0])],
        Float64[],
        Float64[]
    ]

    true_r = [
        Vector{Float64}[[-0.2, -2.4, 0.4], [0.4, 3.2, 0.0]],
        Vector{Float64}[],
        Vector{Float64}[]
    ]

    @test length(nnlist) == length(system)
    @test nnlist.j == true_j
    @test all(nnlist.R .≈ true_R)
    @test all(all(ri .≈ true_ri) for (ri, true_ri) in zip(nnlist.r, true_r))
end
