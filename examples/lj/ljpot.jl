# Find and add source path to Julia search path
cdir = pwd(); ii = findlast("InteratomicPotentials.jl", cdir);
sourcepath = cdir[1:ii[end]] * "/";    
include(sourcepath * "src/setpath.jl");

using Revise, DelimitedFiles, Preprocessing, Potential

app = Preprocessing.initializeapp(cdir,sourcepath);

app.unitstyle = "metal"; # unit system
app.pbc = [1 1 1];       # periodic boundary conditions
app.appname = "lj1";

# data input file 
app.datapath = sourcepath * "examples/data";
app.datafile = "ortho.config"; # configuration file
app.dataformat = "mdp";  # mdp, lammps, json, extxyz
app.datamode = 1;      # 0 => binary, 1 => text   
config = Preprocessing.readconfigfile(app);

# parametric potential descriptors
epsilon = 0.01029849;
sigma = 3.4;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
rcut = 8.5;

potential = Array{Any}(undef, 1)
potential[1] = Potential.addpotential("pair", "nonbonded", [0], "ljpotential", [A; B], rcut);
for i = 1:length(potential)
    include(potential[i].potentialfunction * ".jl");  # include the potential file
end

rcutmax = rcut; eta = [0]; kappa = [0];
etot, e, f = Potential.empiricalpotential(config.x, config.q, config.t, config.a, config.b, config.c, config.pbc, rcutmax, eta, kappa, potential);




