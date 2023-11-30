using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic
using LinearAlgebra


## Settings
num_elements = 1
num_atoms = 100
twojmax = 6
rcutfac = 1.2
radii = [4.0]
rcut0 = 0.989
weight = [1.0]

chem_flag = false
bzero_flag = false
bnorm_flag = false
snap = SNAP(num_atoms, twojmax, [:Ar], rcutfac, 0.00, rcut0, radii, weight, chem_flag, bzero_flag, bnorm_flag)

num_coeffs = length(snap)

print_flag = false


# Create large configuration with num_atoms
min_d = 1.0u"Å"
θ(r) = acos(1.0 - 0.5 * min_d^2 / r^2 )
r = [] 
ri = 1.25u"Å"
n = 0
while n < num_atoms
    ϕ_min = θ(ri)
    ϕ_grid = ϕ_min:2ϕ_min:(π-ϕ_min)
    
    for ϕi in ϕ_grid 
        z = ri * cos(ϕi)
        rj = sqrt(ri^2 - z^2)
        θ_min = θ(rj)
        θ_grid = θ_min:2θ_min:(2*pi - θ_min) 
        for θi in θ_grid 
            x = ri * sin(ϕi) * cos(θi)
            y = ri * sin(ϕi) * sin(θi)
            push!(r, [x, y, z])
            global n+=1
            if n >= num_atoms
                break
            end
        end
        if n >= num_atoms 
            break
        end
    end
    global ri += 1.25u"Å"
end

atoms = [Atom(:Ar, ri) for ri in r]

a = maximum(norm(ri) for ri in r)
box = [[a+1.0*u"Å", 0.0*u"Å", 0.0*u"Å"], [0.0*u"Å", a+1.0*u"Å", 0.0*u"Å"], [0.0*u"Å", 0.0*u"Å", a+1.0*u"Å"]] 
bcs = [DirichletZero(), DirichletZero(), DirichletZero()]

## This is the julia form of the configuration
system = FlexibleSystem(atoms, box, bcs)

# Compute SNAP descriptors
@timed B = compute_local_descriptors(system, snap)
