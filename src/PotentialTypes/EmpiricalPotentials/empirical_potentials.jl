################################################################################
# Types of Potentials
################################################################################

include("lj.jl")
include("bm.jl")
include("coulomb.jl")
include("zbl.jl")

################################################################################
# InteratomicPotentials API implmentations for emperical potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    e = 0.0
    f = fill(SVector{3}(zeros(3)), length(A))
    for ii in 1:length(A)
        for (jj, r, R) in zip(nnlist.j[ii], nnlist.r[ii], nnlist.R[ii])
            e += potential_energy(R, p)
            fo = force(R, r, p)
            f[ii] = f[ii] + fo
            f[jj] = f[jj] - fo
        end
    end
    (; e, f)
end

function potential_energy(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    e = 0.0
    for ii in 1:length(A)
        for R in nnlist.R[ii]
            e += potential_energy(R, p)
        end
    end
    e
end

force(r::SVector{3,<:AbstractFloat}, p::EmpiricalPotential) = force(norm(r), r, p)

function force(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    f = fill(SVector{3}(zeros(3)), length(A))
    for ii in 1:length(A)
        for (jj, r, R) in zip(nnlist.j[ii], nnlist.r[ii], nnlist.R[ii])
            fo = force(R, r, p)
            f[ii] = f[ii] + fo
            f[jj] = f[jj] - fo
        end
    end
    f
end

function virial_stress(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    v = SVector{6}(zeros(6))
    for ii in 1:length(A)
        for (r, R) in zip(nnlist.r[ii], nnlist.R[ii])
            fo = force(R, r, p)
            vi = r * fo'
            v = v + [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]]
        end
    end
    v
end
