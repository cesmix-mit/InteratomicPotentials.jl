# Find and add MDP's path to Julia search path
cdir = pwd(); ii = findlast("MDP", cdir);
if ii===nothing
    display("MDP's path is not found. Please uncomment the above line and set the path to MDP.");
else
    MDPpath = cdir[1:ii[end]] * "/";    
    include(MDPpath * "Installation/setpath.jl");
    push!(LOAD_PATH, pwd());
end

using Revise, DelimitedFiles, Preprocessing, Potential, SnapCpp

app = Preprocessing.initializeapp(cdir,MDPpath);

app.unitstyle = "metal"; # unit system
app.pbc = [1 1 1];       # periodic boundary conditions
app.appname = "lj1";

# data input file 
app.datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MDP/examples/data";
app.datafile = "tiltbox.config"; # configuration file
app.dataformat = "mdp";  # mdp, lammps, json, extxyz
app.datamode = 1;      # 0 => binary, 1 => text   
config = Preprocessing.readconfigfile(app);

# # parametric potential descriptors
# epsilon = 0.01029849;
# sigma = 3.4;
# A = 4*epsilon*sigma^12;
# B = 4*epsilon*sigma^6;
# rcut = 8.5;

# potential = Array{Any}(undef, 1)
# potential[1] = Potential.addpotential("pair", "nonbonded", [0], "ljpotential", [A; B], rcut);
# for i = 1:length(potential)
#     include(potential[i].potentialfunction * ".jl");  # include the potential file
# end

# rcutmax = rcut; eta = [0]; kappa = [0];
# etot, e, f = Potential.empiricalpotential(config.x, config.q, config.t, config.a, 
#                         config.b, config.c, config.pbc, rcutmax, eta, kappa, potential);

libpath = "cpuSnap.dylib"
SnapCpp.loadsnap(libpath)

twojmax = Int32(6)
ntypes = Int32(1)
ncoeffall = Int32(31)
rcutfac = 4.67637;  
rfac0 = 0.99363;
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
snapcoeff = zeros(Float64, ntypes*ncoeffall)

sna = SnapCpp.initsna(snaparam, elemradius, elemweight, snapcoeff)
e, f, v = SnapCpp.snappotential(config.x, config.t, config.a, config.b, config.c, config.pbc, sna)

sna = SnapCpp.initsna(snaparam, elemradius, elemweight)
bi, bd, bv = SnapCpp.snapdescriptors(config.x, config.t, config.a, config.b, config.c, config.pbc, sna)




