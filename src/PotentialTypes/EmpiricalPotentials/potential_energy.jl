function potential_energy(A::AbstractSystem, p::EmpiricalPotential)
    e = 0.0
    nnlist = neighborlist(A, p.rcutoff)

    for ii in 1:length(A)
        for R in nnlist.R[ii]
            e += potential_energy(R, p)
        end
    end
    e
end
