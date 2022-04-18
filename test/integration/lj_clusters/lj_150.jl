@testset "Lennard Jones 150 Cluster Test" begin 
    ϵ = 1.0u"eV"
    σ = 1.0u"Å"
    rcutoff = 100.0u"Å"
    species = [:Ar, :H]
    p = LennardJones(ϵ, σ, rcutoff, species)

    box = [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]] * 1u"Å"
    bcs = [DirichletZero(), DirichletZero(), DirichletZero()]
    l = readlines("integration/lj_clusters/150.xyz")
    r = [ parse.(Float64, split(li)) for li in l ]
    atoms = [Atom(:Ar, ri * u"Å") for ri in r]
    system = FlexibleSystem(atoms, box, bcs)

    e_a_tol = auconvert(1e-5 * 1u"eV")
    es, fs = energy_and_force(system, p)
    @test es ≈ auconvert(-893.310258 * 1u"eV") atol = e_a_tol
    fi = auconvert(0.0 * 1u"eV/Å")
    f_a_tol = auconvert(1e-10 * 1u"eV/Å")
    @test sum(fs)[1] ≈ fi  atol=f_a_tol
    @test sum(fs)[2] ≈ fi  atol=f_a_tol
    @test sum(fs)[3] ≈ fi  atol=f_a_tol

    for ff in fs
        for ffi in ff
            @test ffi ≈ fi atol = f_a_tol
        end
    end
    @test virial(system, p) ≈ auconvert(0.0 * 1u"eV") atol=e_a_tol

end