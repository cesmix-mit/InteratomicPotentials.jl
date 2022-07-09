#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct APPStruct

    sourcepath::String;  # path to source code
    currentdir::String;  # current directory
    codename::String; # MDP
    version::String;  # Set version
    appname::String;  # application name
    platform::String; # CPU or GPU
    cpucompiler::String; # Path to CPU compiler
    mpicompiler::String; # Path to MPI compiler
    gpucompiler::String; # Path to GPU compiler
    mpirun::String;      # Path to MPI run command and MPI run options
    cpuflags::String;    # options for CPU compiler
    gpuflags::String;    # options for GGU compiler
    cpumacros::String;    # options for CPU compiler
    gpumacros::String;    # options for GGU compiler
    cuda_arch::String; 
    enzyme::String; 
    unitstyle::String; 
    ensemblemode::String; 
    descriptor::String; 
    globaloutputs::String; 
    peratomoutputs::String; 
    potentialfile::String;# APP model file name
    configfile::String;# APP model file name    
    weightfile::String;# APP model file name
    model::String;# used to indicate APP model
    buildexec::Int64;
    modelfile::String;# APP model file name    
    modelnumber::Int64;    
    preprocessmode::Int64; # preprocessing mode
    configmode::Int64;
    weightmode::Int64;
    mpiprocs::Int64; # number of MPI ranks
    globalfreq::Int64;
    peratomfreq::Int64;
    unitstylenum::Int64; 
    ensemblemodenum::Int64; 
        
    datapath::String;
    datafile::String;
    dataformat::String;
    datamode::Int64;

    traininglist::Array{Int64,2}; # a list of configurations for training the potential
    validatelist::Array{Int64,2}; # a list of configurations for validating the potential

    dim::Int64; # physical dimension
    nconfigs::Int64; # number of configurations
    maxnumneighbors::Int64; # maximum number of neighbors allowed
    ntimesteps::Int64;# number of time steps
    nab::Int64;# number of atoms per block
    backend::Int64;# backend computing platform
    training::Int64;  # 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    dftdata::Int64;   # 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
    runMD::Int64;        # 0 no MD simulation, 1 -> run MD simulation
    potentialform::Int64;# 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
    #neighpair::Int64;    # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
    #neighcell::Int64;    # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    decomposition::Int64;# 0 -> force decomposition, 1 -> atom decomposition
    energycal::Int64;    # turns energy calculation on or off
    forcecal::Int64;     # turns force calculation on or off
    stresscal::Int64;    # turns stress calculation on or off

    neighskin::Float64;
    neighcell::Int64;       # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    neighmax::Int64;      # maximum number of neighbors allowed
    neighevery::Int64;      # perform neighbor rebuild list for "every" iterations
    neighdelay::Int64;      # delay neighbor rebuild
    neighcheck::Int64;      # 0 -> neighbor rebuild by delay and every settings, 
                             # 1 -> neighbor rebuild occurs if some atom has moved more than half the skin distance
    neighpair::Int64;       # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials

    time::Float64;       # initial time
    dt::Float64;         # time step
    
    coeff::Array{Float64,2};   # coefficients for the potential
    rcutsqmax::Float64;  # square of maximum cutoff radius
    boxoffset::Array{Float64,2}; # offset for simulation box to account for periodic boundary conditions   
    bcs::Array{Int64,2}; # boundary conditions
    pbc::Array{Int64,2};  # periodic boundary conditions
    
    we::Float64;  # energy weight
    wf::Float64;  # force weight

    # atom types
    natomtype::Int64;  # number of atom types
    atomspecies::Array{String,1};
    atommasses::Array{Float64,2};
    atomcharges::Array{Float64,2};
    atomnumbers::Array{Int64,2};
    
    # molecue types
    nmoletype::Int64;  # number of molecule types
    moleinfo::Array{Int64,2}; # information about each molecule type
    
    # simulation and solver parameters
    ndims::Array{Int64,2};     # a list of integers indicating dimensions 
    flag::Array{Int64,2};      # a list of flags 
    simparam::Array{Float64,2}; # simulation parameters
    solparam::Array{Float64,2}; # solver parameters
    nveparam::Array{Float64,2}; # NVE parameters
    nvtparam::Array{Float64,2}; # NVT parameters
    createvelocity::Array{Float64,2}; # NVT parameters
    setatomtypefraction::Array{Float64,2};

    snaprcutfac::Float64;  
    snaprfac0::Float64;  
    snaprmin0::Float64;  
    snaptwojmax::Int64;  
    snapbzeroflag::Int64;  
    snapquadraticflag::Int64;  
    snapswitchflag::Int64;  
    snapchemflag::Int64;  
    snapbnormflag::Int64;  
    snapwselfallflag::Int64;  
    snapnelem::Int64;  
    snapncoeff::Int64;  
    snaparam::Array{Float64,2}; 
    snapelemradius::Array{Float64,2}; 
    snapelemweight::Array{Float64,2}; 
    snapcoeff::Array{Float64,2}; 
    snapcoefffile::String;  
    
    # machine learning potentials
    muml::Array{Float64,2};   # coefficients of ML potential
    K::Int64;          # degree of radial basis functions
    L::Int64;          # degree of spherical harmonics     
    rcutml::Float64;    # cut-off radius for machine learning potential
    descriptornum::Int64; # descriptor flag: 0 -> Spherical Harmonics Bessel
    spectrum::Int64;   # spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    chemtype::Int64;   # 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
    ncmuml::Int64;     # length of mu0
    
    # empirical potential parameters
    eta::Array{Float64,2};   # hyperparameters for emprical potentials
    kappa::Array{Int64,2}; # flag parameters for emprical potentials
    
    # nonbonded single potentials
    np1a::Int64;           # number of nonbonded single potentials 
    pot1a::Array{Int64,2};     # list of nonbonded single potentials 
    mu1a::Array{Float64,2};  # parameters for nonbonded single potentials 
    ncmu1a::Int64;         # length of mu1a
    
    # bonded single potentials
    np1b::Int64;           # number of bonded single potentials 
    pot1b::Array{Int64,2};      # list of bonded single potentials 
    mu1b::Array{Float64,2};   # parameters for bonded single potentials 
    ncmu1b::Int64;         # length of mu1a
    atom1b::Array{Int64,2};     # list of bonded atom types
    
    # nonbonded pair potentials
    np2a::Int64;           # number of nonbonded pair potentials 
    pot2a::Array{Int64,2};     # list of nonbonded pair potentials 
    mu2a::Array{Float64,2};  # parameters for all nonbonded pair potentials 
    ncmu2a::Int64;         # length of mu2a
    rcut2a::Array{Float64,2}; # cut-off radius for each nonbonded pair potential
    
    # bonded pair potentials
    np2b::Int64;               # number of bonded pair potentials 
    pot2b::Array{Int64,2};       # list of bonded pair potentials 
    mu2b::Array{Float64,2};  # parameters for all bonded pair potentials 
    ncmu2b::Int64;             # length of mu2b
    rcut2b::Array{Float64,2};   # cut-off radius for each nonbonded pair potential
    atom2b::Array{Int64,2}; # list of bonded atom pair types
    
    # two-body bond order potentials
    np2c::Int64;               # number of two-body bond order potentials
    pot2c::Array{Int64,2};         # list of two-body bond order potentials
    mu2c::Array{Float64,2};  # parameters for all two-body bond order potentials
    ncmu2c::Int64;             # length of mu2c
    rcut2c::Array{Float64,2};    # cut-off radius for each two-body bond order potential
    atom2c::Array{Int64,2};        # list of bond porder atom types
    
    # nonbonded triplet potentials
    np3a::Int64;           # number of nonbonded triplet potentials
    pot3a::Array{Int64,2};     # list of nonbonded triplet potentials
    mu3a::Array{Float64,2};  # parameters for all nonbonded triplet potentials
    ncmu3a::Int64;         # length of mu3a
    rcut3a::Array{Float64,2}; # cut-off radius for each nonbonded triplet potential
    
    # bonded triplet potentials
    np3b::Int64;               # number of bonded triplet potentials 
    pot3b::Array{Int64,2};         # list of bonded triplet potentials 
    mu3b::Array{Float64,2};  # parameters for all bonded triplet potentials 
    ncmu3b::Int64;             # length of mu3b
    rcut3b::Array{Float64,2};    # cut-off radius for each nonbonded triplet potential
    atom3b::Array{Int64,2}; # list of bonded atom triplet types
    
    # three-body bond order potentials
    np3c::Int64;               # number of three-body bond order potentials
    pot3c::Array{Int64,2};         # list of three-body bond order potentials
    mu3c::Array{Float64,2};  # parameters for all three-body bond order potentials
    ncmu3c::Int64;             # length of mu3c
    rcut3c::Array{Float64,2};    # cut-off radius for each three-body bond order potential
    atom3c::Array{Int64,2};   # list of three-body bond order atom types
    
    # nonbonded quadruplet potentials
    np4a::Int64;           # number of nonbonded quadruplet potentials
    pot4a::Array{Int64,2};     # list of nonbonded quadruplet potentials
    mu4a::Array{Float64,2};  # parameters for all nonbonded quadruplet potentials
    ncmu4a::Int64;          # length of mu4a
    rcut4a::Array{Float64,2}; # cut-off radius for each nonbonded quadruplet potential
    
    # bonded quadruplet potentials
    np4b::Int64;               # number of bonded quadruplet potentials 
    pot4b::Array{Int64,2};         # list of bonded quadruplet potentials 
    mu4b::Array{Float64,2};  # parameters for all bonded quadruplet potentials 
    ncmu4b::Int64;             # length of mu3b
    rcut4b::Array{Float64,2};    # cut-off radius for each nonbonded quadruplet potential
    atom4b::Array{Int64,2}; # list of bonded atom quadruplet types
    
    ncx::Int64; # number of compoments of x
    ncv::Int64; # number of compoments of v
    nce::Int64; # number of compoments of e
    ncf::Int64; # number of compoments of f
    ncq::Int64; # number of compoments of q
    ncmu::Int64; # number of compoments of mu
    nceta::Int64; # number of compoments of eta
    nckappa::Int64; # number of compoments of kappa
            
    lattice;
    domain;
    region;

    APPStruct() = new();
end

function initializeapp(cdir,sourcepath,version="version0")
    app = APPStruct();

    app.sourcepath = sourcepath;
    app.currentdir = cdir;
    app.codename = "MDP";    
    app.version = version;
    app.appname = "app";
    app.platform = "cpu";
    app.cpucompiler = "clang++";
    app.mpicompiler = "mpicxx";
    app.gpucompiler = "clang++";
    app.cuda_arch = "sm_60";
    app.enzyme = "";
    app.mpirun = "mpirun";
    app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
    app.gpuflags = "-lcudart -lcublas";
    app.cpumacros = "";
    app.gpumacros = "";
    app.unitstyle = "";    
    app.descriptor = ""; 
    app.ensemblemode = ""; 
    app.potentialfile = "";
    app.configfile = "";
    app.weightfile = "";   
    app.buildexec = 0;
    app.configmode = 0;    
    app.weightmode = 0;
    app.preprocessmode = 1;
    app.mpiprocs = 1;        

    app.datapath = "";
    app.datafile = "";
    app.dataformat = "";
    app.datamode = 1;

    app.traininglist = reshape([], 0, 2);
    app.validatelist = reshape([], 0, 2);
    
    app.dim = 3;         # physical dimension
    app.nconfigs = 1;    # number of configurations    
    app.maxnumneighbors = 10; # maximum number of neighbors allowed
    app.ntimesteps = 0;  # number of time steps
    app.nab = 4096;      # number of atoms per block
    app.backend = 0;     # backend computing platform
    app.training = 0;     # 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    app.dftdata = 0;      # 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
    app.runMD = 0;        # 0 no MD simulation, 1 -> run MD simulation
    app.potentialform = 0;# 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
    app.decomposition = 0;# 0 -> force decomposition, 1 -> atom decomposition
    app.energycal = 0;    # turns energy calculation on or off
    app.forcecal = 0;     # turns force calculation on or off
    app.stresscal = 0;    # turns stress calculation on or off
    app.time = 0.0;       # initial time
    app.dt = 0.0;         # time step
    app.rcutsqmax = 0.0;  # square of maximum cutoff radius
    app.boxoffset = [0.0 0.0 0.0]; # offset for simulation box to account for periodic boundary conditions
    
    app.globalfreq = 1;       # frequency to print and save default global quantities (pe, ke, ce, te, temp, press)
    app.peratomfreq = 1;     # frequency to save default peratom quantities (id, t, x, f)
    app.globaloutputs = "";    # additional global outputs (stresses, com, vcm)
    app.peratomoutputs = "";   # additional peratom outputs (velocity, mass, virial, image)
    
    app.neighskin = 0.0;
    app.neighcell = 0;       # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    app.neighmax = 100;      # maximum number of neighbors allowed
    app.neighevery = 1;      # perform neighbor rebuild list for "every" iterations
    app.neighdelay = 0;      # delay neighbor rebuild
    app.neighcheck = 1;      # 0 -> neighbor rebuild by delay and every settings, 
                             # 1 -> neighbor rebuild occurs if some atom has moved more than half the skin distance
    app.neighpair = 0;       # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
    
    app.coeff = reshape([], 0, 2);
    app.bcs = [0 0 0 0 0 0]; # boundary conditions
    app.pbc = [1 1 1];       # periodic boundary conditions
    
    app.we = 1.0; # energy weight
    app.wf = 1.0; # force weight

    # atom types
    app.natomtype = 0;  # number of atom types
    app.atommasses = reshape([], 0, 2);
    app.atomcharges = reshape([], 0, 2);
    app.atomnumbers = reshape([], 0, 2);
    
    # molecue types
    app.nmoletype = 0;  # number of molecule types
    app.moleinfo = reshape([], 0, 2); # information about each molecule type
    
    # simulation and solver parameters
    app.ndims = reshape([], 0, 2);     # a list of integers indicating dimensions 
    app.flag = reshape([], 0, 2);      # a list of flags 
    app.simparam = reshape([], 0, 2); # simulation parameters
    app.solparam = reshape([], 0, 2); # solver parameters
    app.nveparam = reshape([], 0, 2); # NVE parameters
    app.nvtparam = reshape([], 0, 2); # NVT parameters
    app.createvelocity = reshape([], 0, 2); # NVT parameters
    app.setatomtypefraction = reshape([], 0, 2);

    # machine learning potentials
    app.muml = reshape([], 0, 2);   # coefficients of ML potential
    app.K = 0;          # degree of radial basis functions
    app.L = 0;          # degree of spherical harmonics     
    app.rcutml = 0.0;    # cut-off radius for machine learning potential
    app.spectrum = 2;   # spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    app.chemtype = 0;   # 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
    app.ncmuml = 0;     # length of mu0
    
    app.snaprcutfac = 0;  
    app.snaptwojmax = 0;
    app.snaprfac0 = 0;
    app.snaprmin0 = 0;
    app.snapbzeroflag = 1;
    app.snapquadraticflag = 0;
    app.snapswitchflag = 1;
    app.snapchemflag = 0;
    app.snapbnormflag = 0;
    app.snapwselfallflag = 0;
    app.snapnelem = 0;    # # of elements
    app.snapncoeff = 0;  # # of coeffients per element
    app.snaparam = reshape([], 0, 2); # snap parameters
    app.snapelemradius = reshape([], 0, 2); # radius per element
    app.snapelemweight = reshape([], 0, 2);   # weight per element
    app.snapcoeff = reshape([], 0, 2);
    app.snapcoefffile = "";    

    # empirical potential parameters
    app.eta = reshape([], 0, 2);   # hyperparameters for emprical potentials
    app.kappa = reshape([], 0, 2); # flag parameters for emprical potentials
    
    # nonbonded single potentials
    app.np1a = 0;           # number of nonbonded single potentials 
    app.pot1a = reshape([], 0, 2);     # list of nonbonded single potentials 
    app.mu1a = reshape([], 0, 2);  # parameters for nonbonded single potentials 
    app.ncmu1a = 0;         # length of mu1a
    
    # bonded single potentials
    app.np1b = 0;           # number of bonded single potentials 
    app.pot1b = reshape([], 0, 2);      # list of bonded single potentials 
    app.mu1b = reshape([], 0, 2);   # parameters for bonded single potentials 
    app.ncmu1b = 0;         # length of mu1a
    app.atom1b = reshape([], 0, 2);     # list of bonded atom types
    
    # nonbonded pair potentials
    app.np2a = 0;           # number of nonbonded pair potentials 
    app.pot2a = reshape([], 0, 2);     # list of nonbonded pair potentials 
    app.mu2a = reshape([], 0, 2);  # parameters for all nonbonded pair potentials 
    app.ncmu2a = 0;         # length of mu2a
    app.rcut2a = reshape([],0,2); # cut-off radius for each nonbonded pair potential
    
    # bonded pair potentials
    app.np2b = 0;               # number of bonded pair potentials 
    app.pot2b = reshape([], 0, 2);       # list of bonded pair potentials 
    app.mu2b = reshape([], 0, 2);  # parameters for all bonded pair potentials 
    app.ncmu2b = 0;             # length of mu2b
    app.rcut2b = reshape([],0,2);   # cut-off radius for each nonbonded pair potential
    app.atom2b = reshape([], 0, 2); # list of bonded atom pair types
    
    # two-body bond order potentials
    app.np2c = 0;               # number of two-body bond order potentials
    app.pot2c = reshape([], 0, 2);         # list of two-body bond order potentials
    app.mu2c = reshape([], 0, 2);  # parameters for all two-body bond order potentials
    app.ncmu2c = 0;             # length of mu2c
    app.rcut2c = reshape([],0,2);    # cut-off radius for each two-body bond order potential
    app.atom2c = reshape([], 0, 2);        # list of bond porder atom types
    
    # nonbonded triplet potentials
    app.np3a = 0;           # number of nonbonded triplet potentials
    app.pot3a = reshape([], 0, 2);     # list of nonbonded triplet potentials
    app.mu3a = reshape([], 0, 2);  # parameters for all nonbonded triplet potentials
    app.ncmu3a = 0;         # length of mu3a
    app.rcut3a = reshape([],0,2); # cut-off radius for each nonbonded triplet potential
    
    # bonded triplet potentials
    app.np3b = 0;               # number of bonded triplet potentials 
    app.pot3b = reshape([], 0, 2);         # list of bonded triplet potentials 
    app.mu3b = reshape([], 0, 2);  # parameters for all bonded triplet potentials 
    app.ncmu3b = 0;             # length of mu3b
    app.rcut3b = reshape([],0,2);    # cut-off radius for each nonbonded triplet potential
    app.atom3b = reshape([], 0, 2); # list of bonded atom triplet types
    
    # three-body bond order potentials
    app.np3c = 0;               # number of three-body bond order potentials
    app.pot3c = reshape([], 0, 2);         # list of three-body bond order potentials
    app.mu3c = reshape([], 0, 2);  # parameters for all three-body bond order potentials
    app.ncmu3c = 0;             # length of mu3c
    app.rcut3c = reshape([],0,2);    # cut-off radius for each three-body bond order potential
    app.atom3c = reshape([], 0, 2);   # list of three-body bond order atom types
    
    # nonbonded quadruplet potentials
    app.np4a = 0;           # number of nonbonded quadruplet potentials
    app.pot4a = reshape([], 0, 2);     # list of nonbonded quadruplet potentials
    app.mu4a = reshape([], 0, 2);  # parameters for all nonbonded quadruplet potentials
    app.ncmu4a = 0;          # length of mu4a
    app.rcut4a = reshape([],0,2); # cut-off radius for each nonbonded quadruplet potential
    
    # bonded quadruplet potentials
    app.np4b = 0;               # number of bonded quadruplet potentials 
    app.pot4b = reshape([], 0, 2);         # list of bonded quadruplet potentials 
    app.mu4b = reshape([], 0, 2);  # parameters for all bonded quadruplet potentials 
    app.ncmu4b = 0;             # length of mu3b
    app.rcut4b = reshape([],0,2);    # cut-off radius for each nonbonded quadruplet potential
    app.atom4b = reshape([], 0, 2); # list of bonded atom quadruplet types
    
    app.ncx = 0; # number of compoments of x
    app.ncv = 0; # number of compoments of v
    app.nce = 0; # number of compoments of e
    app.ncf = 0; # number of compoments of f
    app.ncq = 0; # number of compoments of q
    app.ncmu = 0; # number of compoments of mu
    app.nceta = 0; # number of compoments of eta
    app.nckappa = 0; # number of compoments of kappa
    
    app.lattice = [];
    app.region = [];
    app.domain = [];
    
    return app;

end
