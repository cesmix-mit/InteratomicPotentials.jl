# Contributions by: Stefan Bringuier
################################################################################
# Types of Potentials
################################################################################

include("suttonchen.jl")

export SuttonChen
export get_parameters, get_hyperparameters, serialize_parameters, serialize_hyperparameters

################################################################################
# InteratomicPotentials API implmentations for EAM potentials
################################################################################

function energy_and_force(A::AbstractSystem, p::EmbeddedAtomPotential)
    nnlist = neighborlist(A, p.rcutoff)

    e = 0.0
    f = fill(SVector{3}(zeros(3)), length(A))
    ρ = zeros(length(A))
    dρ = fill(SVector{3}(zeros(3)),length(A))
    for ii in 1:length(A)
        for (jj, r, R) in zip(nnlist.j[ii], nnlist.r[ii], nnlist.R[ii])
            species = unique([atomic_symbol(A, ii), atomic_symbol(A, jj)])
            if (intersect(species, p.species) == species)

                ρ[ii] += rho(R,p)
                ρ[jj] += rho(R,p)
                dρᵢⱼ = drho_dR(R,r,p)
                dρ[ii] = dρ[ii] + dρᵢⱼ
                dρ[jj] = dρ[jj] - dρᵢⱼ

                e += potential_energy_repulsive(R, p)
                fo = force_repulsive(R, r, p)
                f[ii] = f[ii] + fo
                f[jj] = f[jj] - fo
            end
        end
        e -= potential_energy_embedding(ρ[ii],p)
        f[ii] = f[ii] - force_embedding(ρ[ii],dρ[ii],p)
    end
    (; e, f)
end


function force(R::AbstractFloat,r::SVector{3,<:AbstractFloat}, p::EmbeddedAtomPotential)
    ρ,dρ = rho(R,p), drho_dR(R,r,p)
    fo = force_repulsive(R,r,p)
    fo -= force_embedding(ρ,dρ,p)
    return fo
end

force(r::SVector{3,<:AbstractFloat}, p::EmbeddedAtomPotential) = force(norm(r),r,p)

function virial_stress(A::AbstractSystem, p::EmbeddedAtomPotential)
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
