function force(A::AbstractSystem, p::EmpiricalPotential)
    number_of_particles = length(A.particles)
    f = [zeros(3) for i = 1:number_of_particles]
    nnlist = neighborlist(A, p.rcutoff)

    for i in nnlist.i
        j = nnlist.j[i[1]]
        r = nnlist.r[i[1]]
        R = nnlist.R[i[1]]
        for (jj, rj, Rj) in zip(j, r, R)
            fo = force(Rj, rj, p)
            f[i[1]] += fo
            f[jj] -= fo

        end
    end
    return SVector{number_of_particles}([SVector{3}(fi) for fi in f])
end