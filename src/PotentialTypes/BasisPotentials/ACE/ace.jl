using JuLIP: Atoms, JVec, JMat
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
    basis_param :: RPIParams
end

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
    ACE1.forces(get_rpi(rpi), convert_system_to_atoms(A))
end