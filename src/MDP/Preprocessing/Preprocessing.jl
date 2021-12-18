#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Preprocessing

using Revise, JSON

export initializeapp, initializeconfig, initializemdp, preprocessing
export readconfig, readconfigfile, readMDP, readEXTXYZ, readJSON, readJSONfromFolder
export setlattice, setregion, setdomain

include("initializeapp.jl");
include("initializeconfig.jl");
include("initializemdp.jl");
include("readconfig.jl");
include("readconfigfile.jl");
include("readMDP.jl");
include("readEXTXYZ.jl");
include("readJSON.jl");
include("readJSONfromFolder.jl");
include("readweight.jl");
include("writeapp.jl");
include("writeconfig.jl");
include("checkconfig.jl");
include("readLAMMPS.jl");
include("cubemapping.jl");
include("boxperiodicimages.jl");
include("setlattice.jl");
include("setregion.jl");
include("setdomain.jl");
include("writelattice.jl");
include("writeregion.jl");
include("writedomain.jl");

function preprocessing(app)

if app.datafile != "" & app.dataformat != ""
    config = readconfigfile(app);
    app.nconfigs = config.nconfigs;    
elseif app.configmode>0
    # read configurations from data file
    config = readconfig(app, app.configmode);     

    # number of configurations
    app.nconfigs = config.nconfigs;          
else
    app.nconfigs = 0;          
    config = 0;
end

if app.snapcoefffile != ""   
    app.snapcoeff = Main.DelimitedFiles.readdlm(app.snapcoefffile);
end

# nonbonded single potentials
if  length(app.pot1a[:]) > 0
    app.np1a = length(app.pot1a);  # number of nonbonded single potentials 
    app.ncmu1a = length(app.mu1a); # length of mu1a
else
    app.np1a = 0;
    app.ncmu1a = 0;
end

# bonded single potentials
if  length(app.pot1b[:]) > 0
    app.np1b = length(app.pot1b);  # number of nonbonded single potentials 
    app.ncmu1b = length(app.mu1b); # length of mu1a
else
    app.np1b = 0;
    app.ncmu1b = 0;
end

# nonbonded pair potentials
if  length(app.pot2a[:]) > 0
    app.np2a = length(app.pot2a);  # number of nonbonded single potentials 
    app.ncmu2a = length(app.mu2a); # length of mu1a
else
    app.np2a = 0;
    app.ncmu2a = 0;
end

# bonded pair potentials
if  length(app.pot2b[:]) > 0
    app.np2b = length(app.pot2b);  # number of nonbonded single potentials 
    app.ncmu2b = length(app.mu2b); # length of mu1a
else
    app.np2b = 0;
    app.ncmu2b = 0;
end

# two-body bond order potentials
if  length(app.pot2c[:]) > 0
    app.np2c = length(app.pot2c);  # number of nonbonded single potentials 
    app.ncmu2c = length(app.mu2c); # length of mu1a
else
    app.np2c = 0;
    app.ncmu2c = 0;
end

# nonbonded triplet potentials
if  length(app.pot3a[:]) > 0
    app.np3a = length(app.pot3a);  # number of nonbonded single potentials 
    app.ncmu3a = length(app.mu3a); # length of mu1a
else
    app.np3a = 0;
    app.ncmu3a = 0;
end

# bonded triplet potentials
if  length(app.pot3b[:]) > 0
    app.np3b = length(app.pot3b);  # number of nonbonded single potentials 
    app.ncmu3b = length(app.mu3b); # length of mu1a
else
    app.np3b = 0;
    app.ncmu3b = 0;
end

# three-body bond order potentials
if  length(app.pot3c[:]) > 0
    app.np3c = length(app.pot3c);  # number of nonbonded single potentials 
    app.ncmu3c = length(app.mu3c); # length of mu1a
else
    app.np3c = 0;
    app.ncmu3c = 0;
end

# nonbonded quadruplet potentials
if  length(app.pot4a[:]) > 0
    app.np4a = length(app.pot4a);  # number of nonbonded single potentials 
    app.ncmu4a = length(app.mu4a); # length of mu1a
else
    app.np4a = 0;
    app.ncmu4a = 0;
end

# bonded quadruplet potentials
if  length(app.pot4b[:]) > 0
    app.np4b = length(app.pot4b);  # number of nonbonded single potentials 
    app.ncmu4b = length(app.mu4b); # length of mu1a
else
    app.np4b = 0;
    app.ncmu4b = 0;
end

style = app.unitstyle;
if lowercase(style) == "lj" 
    app.unitstylenum = 0;
elseif lowercase(style) == "real"
    app.unitstylenum = 1;    
elseif lowercase(style) == "metal"
    app.unitstylenum = 2;
elseif lowercase(style) == "si"
    app.unitstylenum = 3;
elseif lowercase(style) == "cgs"
    app.unitstylenum = 4;
elseif lowercase(style) == "electron"
    app.unitstylenum = 5;
elseif lowercase(style) == "micro"
    app.unitstylenum = 6;
elseif lowercase(style) == "nano"
    app.unitstylenum = 7;
else
    error("Invalid unit style");
end

# descriptor flag: 0 -> Spherical Harmonics Bessel, 1-> snap
if (lowercase(app.descriptor) == "shp")
    app.descriptornum = 0;
elseif (lowercase(app.descriptor) == "snap") 
    app.descriptornum = 1;    
else
    app.descriptornum = -1; 
end

ensemblemode = app.ensemblemode;
if (lowercase(ensemblemode) == "nve") 
    app.ensemblemodenum = 0;
    app.runMD = 1;  
elseif (lowercase(ensemblemode) == "nvelimit") 
    app.ensemblemodenum = 1;    
    app.runMD = 1;  
elseif (lowercase(ensemblemode) == "nvt") 
    app.ensemblemodenum = 2;        
    app.runMD = 1;  
else
    app.ensemblemodenum = -1;
    app.runMD = 0;  
end

app.natomtype = maximum([length(app.atommasses), length(app.atomspecies)]);
#app.natomtype = length(app.atommasses);
if length(app.atomnumbers)>0
    app.atomnumbers = [0 app.atomnumbers];
else
    app.atomnumbers = reshape([0],1,1);
end
if length(app.atommasses)>0
    app.atommasses = [0.0 app.atommasses];
else
    app.atommasses = reshape([0.0],1,1);
end
if length(app.atomcharges)>0
    app.atomcharges = [0.0 app.atomcharges];
else
    app.atomcharges = reshape([0.0],1,1);
end

app.flag = [app.descriptornum app.spectrum app.training app.runMD app.potentialform app.neighpair app.energycal app.forcecal app.stresscal app.neighcell app.decomposition app.chemtype app.dftdata app.unitstylenum app.ensemblemodenum app.neighcheck];

if length(app.rcut2a) > 0
    rcut2a = app.rcut2a;
else
    rcut2a = reshape([0.0],1,1);
end
if length(app.rcut2b) > 0
    rcut2b = app.rcut2b;
else
    rcut2b = reshape([0.0],1,1);
end
if length(app.rcut2c) > 0
    rcut2c = app.rcut2c;
else
    rcut2c = reshape([0.0],1,1);
end
if length(app.rcut3a) > 0
    rcut3a = app.rcut3a;
else
    rcut3a = reshape([0.0],1,1);
end
if length(app.rcut3b) > 0
    rcut3b = app.rcut3b;
else
    rcut3b = reshape([0.0],1,1);
end
if length(app.rcut3c) > 0
    rcut3c = app.rcut3c;
else
    rcut3c = reshape([0.0],1,1);
end
if length(app.rcut4a) > 0
    rcut4a = app.rcut4a;
else
    rcut4a = reshape([0.0],1,1);
end
if length(app.rcut4b) > 0
    rcut4b = app.rcut4b;
else
    rcut4b = reshape([0.0],1,1);
end
rcut = [reshape([app.rcutml],1,1) rcut2a rcut2b rcut2c rcut3a rcut3b rcut3c rcut4a rcut4b];
rcutmax = maximum(rcut) + app.neighskin;     
app.rcutsqmax = rcutmax^2;        
app.boxoffset = [rcutmax rcutmax rcutmax];
app.simparam = [app.time app.dt];        
app.solparam = [rcutmax app.neighskin];
app.snaparam = [app.snapnelem app.snapncoeff app.snaptwojmax app.snaprcutfac app.snaprfac0 app.snaprmin0 app.snapbzeroflag app.snapswitchflag app.snapquadraticflag app.snapchemflag  app.snapbnormflag app.snapwselfallflag];

app.ndims = zeros(20,1);
app.ndims[1] = app.dim;
app.ndims[2] = app.L;
app.ndims[3] = app.K;
app.ndims[4] = app.ntimesteps;
app.ndims[5] = app.nab;
app.ndims[6] = app.mpiprocs;  
app.ndims[7] = app.backend;
app.ndims[8] = app.nconfigs;
app.ndims[9] = app.natomtype;
app.ndims[10] = app.nmoletype; 
app.ndims[11] = app.neighmax;
app.ndims[12] = app.neighevery;
app.ndims[13] = app.neighdelay;
app.ndims[14] = app.globalfreq;
app.ndims[15] = app.peratomfreq;

if app.configmode>0
    m = 1;
    for i = 1:config.nconfigs
        B2C, C2B = cubemapping(config.a[:,i], config.b[:,i], config.c[:,i]);
        ximages = boxperiodicimages(app.pbc, config.a[:,i], config.b[:,i], config.c[:,i]);    
        n = config.natom[i];
        config.x[:,m:(m+n-1)] = checkconfig(config.x[:,m:(m+n-1)], ximages, B2C, C2B);    
        m = m + n;
    end
end

bindir = "exec";
if app.lattice != []
    writelattice(app.lattice, app.appname * "lattice.bin");
    mv(app.appname * "lattice.bin", app.sourcepath * bindir * "/" * app.appname * "lattice.bin", force = true);
end
if app.region != []
    writeregion(app.region, app.appname * "region.bin");
    mv(app.appname * "region.bin", app.sourcepath * bindir * "/" * app.appname * "region.bin", force = true);    
end
if app.domain != []
    writedomain(app.domain, app.appname * "domain.bin");
    mv(app.appname * "domain.bin", app.sourcepath * bindir * "/" * app.appname * "domain.bin", force = true);
end

if app.nconfigs>0
    if app.training == 0
        config.we = reshape([],0,2);
        config.wf = reshape([],0,2);
    else
        config = readweight(app, config);    
    end    
    writeconfig(config, app.appname * "config.bin");
    mv(app.appname * "config.bin", app.sourcepath * bindir * "/" * app.appname * "config.bin", force = true);
end   

writeapp(app, app.appname * "app.bin");               
mv(app.appname * "app.bin", app.sourcepath * bindir * "/"  * app.appname * "app.bin", force = true);

return app, config

end

end
