using JuLIP: Atoms, JVec, JMat, site_energy, Ref
using JuLIP.Chemistry: AtomicNumber
import ACE1
abstract type ACEPotential <: BasisPotential end 


function convert_system_to_atoms(system::AbstractSystem)
    positions = [ustrip.(position) for position in AtomsBase.position(system)]
    velocities = [ustrip.(v) for v in AtomsBase.velocity(system)]
    masses = ustrip.(AtomsBase.atomic_mass(system))
    atomic_number = AtomicNumber.(AtomsBase.atomic_number(system))
    cell = JMat(ustrip.(vcat(AtomsBase.bounding_box(system)...)))
    pbc = JVec([pbc == Periodic() ? true : false for pbc in AtomsBase.boundary_conditions(system)]...)
    Atoms(X = positions, P = velocities, M = masses, Z = atomic_number, cell = cell, pbc = pbc)
end

struct RPIParams <: BasisParameters
    elements 
    body_order        :: Int
    polynomial_degree :: Int
    wL      
    csp
    r0 
    rcut
end

struct RPI <: ACEPotential
    coefficients 
    basis_params :: RPIParams
end

import Base.length 
length(rpi::RPIParams) = length(get_rpi(rpi))

function get_rpi(rpi::RPIParams)
    ACE1.rpi_basis(;
    species = rpi.elements,
    N       = rpi.body_order - 1,
    maxdeg  = rpi.polynomial_degree,
    D       = ACE1.SparsePSHDegree(; wL = rpi.wL, csp = rpi.csp),
    r0      = rpi.r0, 
    rin     = 0.65*rpi.r0,
    rcut    = rpi.rcut,
    pin     = 0,
    )
end

function evaluate_basis(A::AbstractSystem, rpi::RPIParams)
    ACE1.energy(get_rpi(rpi), convert_system_to_atoms(A))
end

function evaluate_basis_d(A::AbstractSystem, rpi::RPIParams)
    ftemp = ACE1.forces(get_rpi(rpi), convert_system_to_atoms(A))
    f = [zeros(3, length(rpi)) for i = 1:length(A)]

    for i = 1:length(A)
        for j = 1:3 
            for k = 1:length(rpi)
                f[i][j, k] = ftemp[k][i][j]
            end
        end
    end
    f
end

function evaluate_basis_v(A::AbstractSystem, rpi::RPIParams)
    Wtemp = ACE1.virial(get_rpi(rpi), convert_system_to_atoms(A))
    W = zeros(6, length(rpi))

    for k = 1:length(rpi)
        count = 1
        for (i, j) in zip( [1, 2, 3, 3, 3, 2], [1, 2, 3, 2, 1, 1])
            W[count, k] = Wtemp[k][i, j]
            count +=1
        end
    end
    W
end

function evaluate_full(A::AbstractSystem, rpi::RPIParams)
    B = hcat([site_energy(get_rpi(rpi), convert_system_to_atoms(A), i) for i = 1:length(A)]...)'
    dB = reshape(hcat(evaluate_basis_d(A, rpi)...), :, length(rpi))
    W = evaluate_basis_v(A, rpi)
    return B, dB, W
end