@testset "Neighbor List Unit Tests" begin
    position = (@SVector [1.0, 0.0, 0.0])u"bohr"
    atoms = [
        :Ar => 0.5 * position,
        :Ar => 0.75 * position,
        :H => 1.1 * position
    ]
    box = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]u"bohr"
    system = periodic_system(atoms, box)

    rcutoff = 2.0
    nnlist = neighborlist(system, rcutoff)

    nnlist_r = nnlist.r
    true_r = [
        [SVector{3}([0.25, 0.0, 0.0]), SVector{3}([-0.40, 0.0, 0.0])],
        [SVector{3}([0.35, 0.0, 0.0])],
        []
    ]

    for (ri, ri_true) in zip(nnlist_r, true_r)
        for (rij, rij_true) in zip(ri, ri_true)
            @test isapprox(sum(rij + rij_true), 0.0, atol=eps(2.0))
        end
    end
end
