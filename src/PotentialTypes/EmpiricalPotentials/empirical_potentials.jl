################################################################################
# Types of Potentials
################################################################################

include("lj.jl")
include("bm.jl")
include("coulomb.jl")
include("zbl.jl")
include("morse.jl")

export LennardJones, BornMayer, Coulomb, ZBL, Morse

################################################################################
# InteratomicPotentials API implmentations for emperical potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    e = 0.0
    f = fill(SVector{3}(zeros(3)), length(A))
    for ii in 1:length(A)
        for (jj, r, R) in zip(nnlist.j[ii], nnlist.r[ii], nnlist.R[ii])
            species = unique([atomic_symbol(A, ii), atomic_symbol(A, jj)])
            if (intersect(species, p.species) == species)
                e += potential_energy(R, p)
                fo = force(R, r, p)
                f[ii] = f[ii] + fo
                f[jj] = f[jj] - fo
            end
        end
    end
    (; e=e * u"hartree", f=f * u"hartree / bohr")
end

force(r::SVector{3}, p::EmpiricalPotential) = force(norm(r), r, p)

function virial_stress(A::AbstractSystem, p::EmpiricalPotential)
    nnlist = neighborlist(A, p.rcutoff)

    v = @SVector zeros(6)
    for ii in 1:length(A)
        for (r, R) in zip(nnlist.r[ii], nnlist.R[ii])
            vi = r * force(R, r, p)'
            v += @SVector [vi[1, 1], vi[2, 2], vi[3, 3], vi[3, 2], vi[3, 1], vi[2, 1]]
        end
    end
    v * u"hartree"
end
