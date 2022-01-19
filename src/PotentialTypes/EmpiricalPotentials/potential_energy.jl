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
    e
end

function potential_energy2(A::AbstractSystem, p::EmpiricalPotential)
    e = 0.0
    nnlist = neighborlist(A, p.rcutoff)

    for ii in 1:length(A)
        for R in nnlist.R[ii]
            e += potential_energy(R, p)
        end
    end
    e
end
