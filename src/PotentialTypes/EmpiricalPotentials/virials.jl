#############################################################################
############################## Simple Virial ################################
#############################################################################
function virial(A::AbstractSystem, p::EmpiricalPotential)
    v = virial_stress(A, p)
    return v[1] + v[2] + v[3]
end



function virial_stress(A::AbstractSystem, p::EmpiricalPotential)
    v = zeros(6)
    nnlist = neighborlist(A, p.rcutoff)

    for i in nnlist.i
        j = nnlist.j[i[1]]
        r = nnlist.r[i[1]]
        R = nnlist.R[i[1]]
        for (jj, rj, Rj) in zip(j, r, R)
            fo = force(Rj, rj, p)
            vi = rj * fo'
            v += [vi[1,1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]]
        end
    end
    return SVector{6}(v)
end