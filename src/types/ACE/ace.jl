using JuLIP: Atoms, JVec, JMat, site_energy, Ref
using JuLIP.Chemistry: AtomicNumber
import ACE1

include("ace_params.jl")

struct ACE <: BasisPotential 
    coefficients :: AbstractVector
    basis_params :: ACEParams
end 


get_rcutoff(ace::ACE) = ace.basis_params.rcutoff
get_species(ace::ACE) = ace.basis_params.species

function convert_system_to_atoms(system::AbstractSystem)
    positions = [ustrip.(position) for position in AtomsBase.position(system)]
    velocities = [ustrip.(v) for v in AtomsBase.velocity(system)]
    masses = ustrip.(AtomsBase.atomic_mass(system))
    atomic_number = AtomicNumber.(AtomsBase.atomic_number(system))
    cell = JMat(ustrip.(vcat(AtomsBase.bounding_box(system)...)))
    pbc = JVec([pbc == Periodic() ? true : false for pbc in AtomsBase.boundary_conditions(system)]...)
    Atoms(X = positions, P = velocities, M = masses, Z = atomic_number, cell = cell, pbc = pbc)
end

function evaluate_basis(A::AbstractSystem, ace_params::ACEParams)
    ACE1.energy(get_rpi(ace_params), convert_system_to_atoms(A))
end

function evaluate_basis_d(A::AbstractSystem, ace_params::ACEParams)
    ftemp = ACE1.forces(get_rpi(ace_params), convert_system_to_atoms(A))
    f = [zeros(3, length(ace_params)) for i = 1:length(A)]

    for i = 1:length(A)
        for j = 1:3 
            for k = 1:length(ace_params)
                f[i][j, k] = ftemp[k][i][j]
            end
        end
    end
    f
end

function evaluate_basis_v(A::AbstractSystem, ace_params::ACEParams)
    Wtemp = ACE1.virial(get_rpi(ace_params), convert_system_to_atoms(A))
    W = zeros(6, length(ace_params))

    for k = 1:length(ace_params)
        count = 1
        for (i, j) in zip( [1, 2, 3, 3, 3, 2], [1, 2, 3, 2, 1, 1])
            W[count, k] = Wtemp[k][i, j]
            count +=1
        end
    end
    W
end

function evaluate_full(A::AbstractSystem, ace_params::ACEParams)
    B = hcat([site_energy(get_rpi(ace_params), convert_system_to_atoms(A), i) for i = 1:length(A)]...)'
    dB = reshape(hcat(evaluate_basis_d(A, ace_params)...), :, length(ace_params))
    W = evaluate_basis_v(A, ace_params)
    return B, dB, W
end