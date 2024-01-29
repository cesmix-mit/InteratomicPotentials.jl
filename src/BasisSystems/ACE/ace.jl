using JuLIP: Atoms, JVec, JMat, site_energy, Ref
using JuLIP.Chemistry: AtomicNumber
import ACE1
using Unitful, UnitfulAtomic

struct ACE <: BasisSystem
    species           :: Vector{Symbol}
    body_order        :: Int   
    polynomial_degree :: Int
    rcutoff           :: Real
    wL                :: Real
    csp               :: Real
    r0                :: Real 
    rpib              :: ACE1.RPIBasis
end 

function ACE(species, body_order, polynomial_degree, rcutoff;
             wL = 1.5, csp = 1.0, r0 = 2.5)
    rpib = ACE1.rpi_basis(;
                species = species,
                N       = body_order - 1,
                maxdeg  = polynomial_degree,
                D       = ACE1.SparsePSHDegree(; wL = wL, csp = csp),
                r0      = r0, 
                rin     = 0.65*r0,
                rcut    = rcutoff,
                pin     = 0,
           )
    return ACE(species, body_order, polynomial_degree,
               rcutoff, wL, csp, r0, rpib)
end

function ACE(; species = [:X], body_order = 4, polynomial_degree = 6, 
               rcutoff = 5.0, wL = 1.5, csp = 1.0, r0 = 2.5)
    rpib = ACE1.rpi_basis(;
                species = species,
                N       = body_order - 1,
                maxdeg  = polynomial_degree,
                D       = ACE1.SparsePSHDegree(; wL = wL, csp = csp),
                r0      = r0, 
                rin     = 0.65*r0,
                rcut    = rcutoff,
                pin     = 0,
           )
    return ACE(species, body_order, polynomial_degree,
               rcutoff, wL, csp, r0, rpib)
end

function get_rpi(ace)
    return ace.rpib
end

function get_rcutoff(ace::ACE)
    return ace.rcutoff
end

function get_species(ace::ACE)
    return ace.species
end

Base.length(ace::ACE) = length(ace.rpib)

function convert_system_to_atoms(system::AbstractSystem)
    positions = [ustrip.(position) for position in AtomsBase.position(system)]
    velocities = [ustrip.(v) for v in AtomsBase.velocity(system)]
    masses = ustrip.(AtomsBase.atomic_mass(system))
    atomic_number = AtomicNumber.(AtomsBase.atomic_number(system))
    cell = JMat(ustrip.(transpose(hcat(AtomsBase.bounding_box(system)...))))
    pbc = JVec([pbc == Periodic() ? true : false for pbc in AtomsBase.boundary_conditions(system)]...)
    return Atoms(X = positions, P = velocities, M = masses, Z = atomic_number, cell = cell, pbc = pbc)
end

function compute_local_descriptors(A::AbstractSystem, ace::ACE)
    return [ site_energy(ace.rpib, convert_system_to_atoms(A), i)
             for i = 1:length(A)]
end

function compute_force_descriptors(A::AbstractSystem, ace::ACE)
    ftemp = ACE1.forces(ace.rpib, convert_system_to_atoms(A))
    # TODO: settle on efficient and consistent descriptor layout
    f = [ [zeros(length(ace)) for _ in 1:3] for i in 1:length(A)]
    for i = 1:length(A)
        for j = 1:3 
            for k = 1:length(ace)
                f[i][j][k] = ftemp[k][i][j]
            end
        end
    end
    return f
end

function compute_virial_descriptors(A::AbstractSystem, ace::ACE)
    Wtemp = ACE1.virial(ace.rpib, convert_system_to_atoms(A))
    W = zeros(6, length(ace))
    for k = 1:length(ace)
        count = 1
        for (i, j) in zip( [1, 2, 3, 3, 3, 2], [1, 2, 3, 2, 1, 1])
            W[count, k] = Wtemp[k][i, j]
            count +=1
        end
    end
    return W
end

function compute_all_descriptors(A::AbstractSystem, ace::ACE)
    return compute_local_descriptors(A, ace),
           compute_force_descriptors(A, ace),
           compute_virial_descriptors(A, ace)
end
