#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct CONFIGStruct

    dim::Int64;
    ncx::Int64;
    ncv::Int64;
    ncq::Int64;
    nce::Int64;
    ncf::Int64;
    ncs::Int64;
    nci::Int64;
    nct::Int64;
    ncg::Int64;
    ncz::Int64;
    ncm::Int64;
    nco::Int64;
    ncl::Int64;
    ncp::Int64;
    nconfigs::Int64;
    natom::Array{Int64,2};   # number of atoms per each configuration
    natomall::Int64;# number of atoms for all configurations
    
    # simulation box for each configuration
    a::Array{Float64,2}; # the 1st principal vector of the simulation box
    b::Array{Float64,2}; # the 2nd principal vector of the simulation box
    c::Array{Float64,2}; # the 3rd principal vector of the simulation box
    lattice::Array{Float64,2};  # energies for all configurations
    pbc::Array{Int64,2};  # energies for all configurations
    e::Array{Float64,2};  # energies for all configurations
    stress::Array{Float64,2};  # energies for all configurations

    we::Array{Float64,2};  # energies weights for all configurations
    wf::Array{Float64,2};  # energies weights for all configurations
    ws::Array{Float64,2};  # energies weights for all configurations

    t::Array{Int64,2};   # atom types for all configurations
    Z::Array{Int64,2};   # atom types for all configurations
    move::Array{Int64,2};   # atom types for all configurations
    tags::Array{Int64,2};   # atom types for all configurations
    group::Array{Int64,2};   # atom types for all configurations
    mass::Array{Float64,2}; # atom positions for all configurations
    x::Array{Float64,2}; # atom positions for all configurations
    q::Array{Float64,2}; # atom charges for all configurations
    v::Array{Float64,2}; # atom velocities for all configurations
    f::Array{Float64,2}; # atom forces for all configurations    

    CONFIGStruct() = new();
end

function initializeconfig(app)

    config = CONFIGStruct();     
    nconfigs = 0;            # number of configurations
    dim = 0
    nci = 0;
    nct = 0;
    ncg = 0;
    ncx = 0;
    ncv = 0;
    ncf = 0;
    ncq = 0;
    nce = 0;
    ncs = 0;
    ncz = 0;
    ncm = 0;
    nco = 0;
    ncl = 0;
    ncp = 0;

    config.nconfigs = nconfigs;
    config.dim = dim;
    config.ncx = ncx;
    config.ncv = ncv;
    config.ncq = ncq;
    config.nce = nce;
    config.ncf = ncf;
    config.ncs = ncs;
    config.nci = nci;
    config.nct = nct;
    config.ncg = ncg;
    config.ncz = ncz;
    config.ncm = ncm;
    config.nco = nco;
    config.ncl = ncl;
    config.ncp = ncp;
    
    config.natom = ones(1, nconfigs);   # number of atoms per each configuration
    config.natomall = sum(config.natom);# number of atoms for all configurations

    # simulation box for each configuration
    config.a = zeros(dim, nconfigs); # the 1st principal vector of the simulation box
    config.b = zeros(dim, nconfigs); # the 2nd principal vector of the simulation box
    config.c = zeros(dim, nconfigs); # the 3rd principal vector of the simulation box
    config.lattice = zeros(ncl, nconfigs); # stresses for all configurations
    config.pbc = zeros(ncp, nconfigs); # periodic boundary conditions
    config.e = zeros(nce, nconfigs);  # energies for all configurations
    config.stress = zeros(ncs, nconfigs); # stresses for all configurations

    config.we = ones(1, nconfigs);      # energy weight per each configuration
    config.wf = ones(1, nconfigs);      # force weight per each configuration
    config.ws = ones(1, nconfigs);      # stress weight per each configuration

    config.Z = zeros(ncz, config.natomall);    # atom numbers for all configurations
    config.mass = zeros(ncm, config.natomall); # atom masses for all configurations
    config.move = zeros(nco, config.natomall);   # atom move masks for all configurations
    config.tags = zeros(nci, config.natomall); # atom tags for all configurations
    config.t = zeros(nct, config.natomall);   # atom types for all configurations
    config.group = zeros(ncg, config.natomall);   # atom groups for all configurations
    config.x = zeros(ncx, config.natomall); # atom positions for all configurations
    config.q = zeros(ncq, config.natomall); # atom charges for all configurations
    config.v = zeros(ncv, config.natomall); # atom velocities for all configurations
    config.f = zeros(ncf, config.natomall); # atom forces for all configurations
    
    return config;
end
