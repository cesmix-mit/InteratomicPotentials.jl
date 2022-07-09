cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "setpath.jl");

# Read xyz file to obtain configurations 
filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster/curated_lj_cluster.xyz";
species = ["Ar"];
config = Preprocessing.readEXTXYZ(filename, species);

# shift atom positions to put them in the box [0,10]x[0,10]x[0,10]
config.x = config.x .+ 5.0

# cumalative sum of numbers of atoms 
natom = [0; cumsum(config.natom[:])];

ci = 1; # configuration index in the data 
x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
t = config.t[:,ci]                        # atom types of configuration ci
q = []

# 3 vectors define the simulation box
a = [10.0; 0.0; 0.0]
b = [0.0; 10.0; 0.0]
c = [0.0; 0.0; 10.0]

# periodic boundary conditions
pbc = [0; 0; 0]

# LJ parameters
epsilon = 1.0; 
sigma = 1.0;
A = 2*epsilon*sigma^12;
B = 2*epsilon*sigma^6;
rcut = 4.0;

# define a list of potentials
potential = Array{Any}(undef, 1)
potential[1] = Potential.addpotential("pair", "nonbonded", [0], "ljpotential", [A; B], rcut);
for i = 1:length(potential)
    include(potential[i].potentialfunction * ".jl");  # include the potential file
end

# compute total energy, peratom energies, and forces
rcutmax = rcut; eta = [0]; kappa = [0];
etot, e, f = Potential.empiricalpotential(config.x[:,1:13], q, t, a, b, c, pbc, rcutmax, eta, kappa, potential);

# snap parameters
twojmax = Int32(6)
ntypes = Int32(1)
ncoeffall = Int32(31)
rcutfac = 4.0;  
rfac0 = 1.0;
rmin0 = 0;
bzeroflag = 0;
quadraticflag = 0;
switchflag = 1;
chemflag = 0;
bnormflag = 0;
wselfallflag = 0;

snaparam = [ntypes, twojmax, rcutfac, rfac0, rmin0, bzeroflag, switchflag, quadraticflag, chemflag, bnormflag, wselfallflag]
elemradius = 0.5*ones(Float64, ntypes)
elemweight = ones(Float64, ntypes)

# intialize sna struct to store snap parameters 
sna = Potential.initsna(snaparam, elemradius, elemweight)

# bispectrum, bispectrum derivatives, bispectrum virial
bi, bd, bv = Potential.snapdescriptors(x, t, a, b, c, pbc, sna)

# using LAMMPS.jl to compute energy and forces
unitstyle = "lj"
pair_style = "lj/cut $rcut"  
pair_coeff = ["* * $epsilon $sigma"]
lbi, lbd, lbv, re, rf, rv = Potential.lammpssnapdescriptors(x, t, a, b, c, pbc, 
            unitstyle, pair_style, pair_coeff, snaparam, elemradius, elemweight)

# compare with LAMMPS.jl
print("Maximum error in bispectrum components: "); display(maximum(abs.(bi[:] .- lbi[:])));
print("Maximum error in bispectrum derivatives: "); display(maximum(abs.(bd[:] .- lbd[:])));
print("Maximum error in bispectrum virials: "); display(maximum(abs.(bv[:] .- lbv[:])));
print("Maximum error in reference energy: "); display(abs(etot - re));
print("Maximum error in reference forces: "); display(maximum(abs.(f[:] .- rf[:])));




