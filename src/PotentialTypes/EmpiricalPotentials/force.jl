function force(A::AbstractSystem, p::EmpiricalPotential)
    f = fill(zeros(3), length(A))
    nnlist = neighborlist(A, p.rcutoff)

    for ii in 1:length(A)
        j = nnlist.j[ii]
        r = nnlist.r[ii]
        R = nnlist.R[ii]
        for (jj, rj, Rj) in zip(j, r, R)
            fo = force(Rj, rj, p)
            f[ii] += fo
            f[jj] -= fo
        end
    end
    [SVector{3}(fi) for fi in f]
end
