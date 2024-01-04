function load_configuration(filename::String)
    box = [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]] * 1u"Å"
    bcs = [DirichletZero(), DirichletZero(), DirichletZero()]
    l = readlines(filename)
    r = [parse.(Float64, split(li)) for li in l]
    atoms = [Atom(:Ar, ri * u"Å") for ri in r]
    FlexibleSystem(atoms, box, bcs)
end

function generate_lj_cluster(N::Integer)
    min_d = 1.0u"Å"
    θ(r) = acos(1.0 - 0.5 * min_d^2 / r^2)
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
            θ_grid = θ_min:2θ_min:(2π-θ_min)
            for θi in θ_grid
                x = ri * sin(ϕi) * cos(θi)
                y = ri * sin(ϕi) * sin(θi)
                push!(r, [x, y, z])
                n += 1
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

    a = maximum(norm(ri) for ri in r) + 1.0u"Å"
    z = 0.0u"Å"
    box = [[a, z, z], [z, a, z], [z, z, a]]
    bcs = [DirichletZero(), DirichletZero(), DirichletZero()]
    FlexibleSystem(atoms, box, bcs)
end
