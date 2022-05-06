@testset "Lennard Jones 10000 Cluster Test" begin 
    ϵ = 1.0u"eV"
    σ = 1.0u"Å"
    rcutoff = 100.0u"Å"
    species = [:Ar, :H]
    p = LennardJones(ϵ, σ, rcutoff, species)

    min_d = 1.0u"Å"
    θ(r) = acos(1.0 - 0.5 * min_d^2 / r^2 )
    N = 10000
    r = [] 
    ri = 1.25u"Å"
    n = 0
    while n < N
        ϕ_min = θ(ri)
        ϕ_grid = ϕ_min:2ϕ_min:(π-ϕ_min)
        
        for ϕi in ϕ_grid 
            z = ri * cos(ϕi)
            rj = sqrt(ri^2 - z^2)
            θ_min = θ(rj)
            θ_grid = θ_min:2θ_min:(2*pi - θ_min) 
            for θi in θ_grid 
                x = ri * sin(ϕi) * cos(θi)
                y = ri * sin(ϕi) * sin(θi)
                push!(r, [x, y, z])
                n+=1
                if n > N
                    break
                end
            end
            if n > N 
                break
            end
        end
        ri += 1.25u"Å"
    end

    atoms = [Atom(:Ar, ri) for ri in r]

    a = maximum(norm(ri) for ri in r)
    box = [[a+1.0*u"Å", 0.0*u"Å", 0.0*u"Å"], [0.0*u"Å", a+1.0*u"Å", 0.0*u"Å"], [0.0*u"Å", 0.0*u"Å", a+1.0*u"Å"]] 
    bcs = [DirichletZero(), DirichletZero(), DirichletZero()]
    system = FlexibleSystem(atoms, box, bcs)

    @time es, fs = energy_and_force(system, p)
    println("Energy: $(es)")
    @test ustrip(es) isa Float64
    @test es > -1e6 * u"hartree" 
    @test length(fs) == n 
    @test sum(sum(fs)) ≈ auconvert(0.0 * 1u"eV/Å") atol =auconvert(1e-5*u"eV/Å")

end