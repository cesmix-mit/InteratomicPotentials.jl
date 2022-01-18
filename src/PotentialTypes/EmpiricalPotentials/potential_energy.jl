function potential_energy(A::AbstractSystem, p::EmpiricalPotential)
    number_of_particles = length(A.particles)
    e = 0.0
    nnlist = neighborlist(A, p.rcutoff)

    for i in nnlist.i
        R = nnlist.R[i[1]]
        for Ri in R
            e += potential_energy(Ri, p)
        end
    end
    return e
end
