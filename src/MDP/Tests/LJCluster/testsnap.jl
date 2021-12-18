cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "setpath.jl");

using Preprocessing, Potential, SnapCpp 

# load Snap library 
libpath = "/Users/ngoccuongnguyen/GitHub/InteratomicPotentials.jl/src/MDP/SnapCpp/cpuSnap.dylib"
SnapCpp.loadsnap(libpath)

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

# 3 vectors define the simulation box
a = [10.0; 0.0; 0.0]
b = [0.0; 10.0; 0.0]
c = [0.0; 0.0; 10.0]

# periodic boundary conditions
pbc = [0; 0; 0]

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
sna = SnapCpp.initsna(snaparam, elemradius, elemweight)

# bispectrum, bispectrum derivatives, bispectrum virial
bi, bd, bv = SnapCpp.snapdescriptors(x, t, a, b, c, pbc, sna)


