# Contributions by: Stefan Bringuier

using AtomsBase
using InteratomicPotentials
using StaticArrays
using Unitful
using UnitfulAtomic

# Reference: A. P. Sutton and J. Chen, Phil. Mag. Lett. 61, 139 (1990)
ϵ = 1.2832e-2
a,c = 3.61,39.432
nn,m = 9,6
rcutoff = 10.0
rc = 5.5
sc = SuttonChen(ϵ,a,c,nn,m,rc, [:Cu,:Cu])

pair_dist,energy = [], []
lcell = 15.000
for i in 0.0:0.05:5.0
    atom1 = Atom(element, (@SVector [0.0, 0.0, 0.0]) * u"Å")
    atom2 = Atom(element, (@SVector [1.0+i, 0.0, 0.0]) * u"Å")
    atoms = [atom1, atom2]
    box = [[lcell, 0.0, 0.0], [0.0, lcell, 0.0], [0.0, 0.0, lcell]]
    bcs = [Periodic(), Periodic(), Periodic()]
    system = FlexibleSystem(atoms, box * u"Å", bcs)
    n = length(system)
    nnlist = neighborlist(system, rcutoff)
    push!(pair_dist,nnlist.R[1][1])
    push!(energy,potential_energy(system,sc))
end

min_energy, iloc = findmin(energy/length(system))
min_dist = pair_dist[iloc]
@test -1.1 < min_energy < -1.3
@test 2.05e0 < min_dist < 2.2e0

e_binding = -3.54 # eV/atom
eqpair_dist = 2.578 # Å
lat_param = 3.615 # Å

function fcc_cu(;latp=3.615)
    atom1 = Atom(element, (@SVector [0.0, 0.0, 0.0]) * latp*u"Å")
    atom2 = Atom(element, (@SVector [0.5, 0.5, 0.0]) * latp*u"Å")
    atom3 = Atom(element, (@SVector [0.5, 0.0, 0.5]) * latp*u"Å")
    atom4 = Atom(element, (@SVector [0.0, 0.5, 0.5]) * latp*u"Å")
    atoms = [atom1,atom2,atom3,atom4]
    box = [[latp, 0.0, 0.0], [0.0, latp, 0.0], [0.0, 0.0, latp]]
    bcs = [Periodic(), Periodic(), Periodic()]
    system = FlexibleSystem(atoms, box * u"Å", bcs)
    return system
end

pair_dist,fcc_energy  = [], []
for latp in 2.25:0.005:8.00
    system = fcc_cu(latp=latp)
    n = length(system)
    nnlist = neighborlist(system, rcutoff)
    push!(pair_dist,nnlist.R[1][1])
    push!(fcc_energy,potential_energy(system,sc))
end

min_fcc_energy, iloc = findmin(fcc_energy/length(system))
min_dist = pair_dist[iloc]
@test -3.545 < min_fcc_energy < -3.550
@test 2.54e0 < min_dist < 2.58e0


force_equil_fcc  = force(fcc_cu(),sc)
@test sum(force_equil_fcc) ≈ SVector(0.0,0.0,0.0)
