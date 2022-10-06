using JuLIP: Atoms, JVec, JMat, site_energy, Ref
using JuLIP.Chemistry: AtomicNumber
import ACE1

struct ACE <: BasisSystem
    species           :: Vector{Symbol}
    body_order        :: Int   
    polynomial_degree :: Int
    wL                :: Real
    csp               :: Real
    r0                :: Real 
    rcutoff           :: Real
end 

function ACE(species :: Vector{Symbol}; body_order = 4, polynomial_degree = 6, 
            wL = 1.0, csp = 1.0, r0 = 0.95, rcutoff = 5.0)
    ACE(species, body_order, polynomial_degree, wL, csp, r0, rcutoff)
end

length(ace::ACE) = length(get_rpi(ace))

function get_rpi(ace::ACE)
    ACE1.rpi_basis(;
    species = ace.species,
    N       = ace.body_order - 1,
    maxdeg  = ace.polynomial_degree,
    D       = ACE1.SparsePSHDegree(; wL = ace.wL, csp = ace.csp),
    r0      = ace.r0, 
    rin     = 0.65*ace.r0,
    rcut    = ace.rcutoff,
    pin     = 0,
    )
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

function compute_local_descriptors(A::AbstractSystem, ace::ACE)
    [site_energy(get_rpi(ace), convert_system_to_atoms(A), i) for i = 1:length(A)]
end

function compute_force_descriptors(A::AbstractSystem, ace::ACE)
    ftemp = ACE1.forces(get_rpi(ace), convert_system_to_atoms(A))
    f = [zeros(3, length(ace)) for i = 1:length(A)]

    for i = 1:length(A)
        for j = 1:3 
            for k = 1:length(ace)
                f[i][j, k] = ftemp[k][i][j]
            end
        end
    end
    f
end

function compute_virial_descriptors(A::AbstractSystem, ace::ACE)
    Wtemp = ACE1.virial(get_rpi(ace), convert_system_to_atoms(A))
    W = zeros(6, length(ace))

    for k = 1:length(ace)
        count = 1
        for (i, j) in zip( [1, 2, 3, 3, 3, 2], [1, 2, 3, 2, 1, 1])
            W[count, k] = Wtemp[k][i, j]
            count +=1
        end
    end
    W
end

function compute_all_descriptors(A::AbstractSystem, ace::ACE)
    compute_local_descriptors(A, ace), compute_force_descriptors(A, ace), compute_virial_descriptors(A, ace)
end