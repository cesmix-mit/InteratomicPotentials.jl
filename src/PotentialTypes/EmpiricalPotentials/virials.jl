function virial(A::AbstractSystem, p::EmpiricalPotential)
    v = virial_stress(A, p)
    v[1] + v[2] + v[3]
end

function virial_stress(A::AbstractSystem, p::EmpiricalPotential)
    v = zeros(6)
    nnlist = neighborlist(A, p.rcutoff)

    for ii in 1:length(A)
        r = nnlist.r[ii]
        R = nnlist.R[ii]
        for (rj, Rj) in zip(r, R)
            fo = force(Rj, rj, p)
            vi = rj * fo'
            v += [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]]
        end
    end
    SVector{6}(v)
end
