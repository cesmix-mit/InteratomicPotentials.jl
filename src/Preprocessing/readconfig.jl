#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function readconfig(app, mode)

filename = app.configfile;

if mode==0 # binary    
    tmp = reinterpret(Float64,read(filename));
else # text    
    tmp = Main.DelimitedFiles.readdlm(filename);
end
#tmp = tmp[:];

if mode==4
    config = readLAMMPS(tmp);
    config.nconfigs = length(config.natom); # number of configurations
    config.natomall = sum(config.natom);# number of atoms for all configurations
    return config;
end

if mode<=2
    nconfigs = int64(tmp[1]);          # number of configurations
    dim = int64(tmp[2]);
    ncq = int64(tmp[3]);
    ncv = int64(tmp[4]);
    nce = int64(tmp[5]);
    ncf = int64(tmp[6]);
    ncx = dim;

    app.nconfigs = nconfigs; # number of configurations
    app.dim = dim; # physical dimension
    app.ncx = dim; # number of compoments of x
    app.ncv = ncv; # number of compoments of v
    app.ncf = ncf; # number of compoments of f
    app.ncq = ncq; # number of compoments of q
    app.nce = nce; # number of compoments of e
    config = initializeconfig(app);
end

if mode<=1
    n2 = 6;
    m  = (1+nce+dim*dim)*nconfigs;
    n1 = n2+1;
    n2 = n2+m;    
    tm = reshape(tmp[n1:n2], nconfigs, (1+nce+dim*dim)); 
    config.natom = tm[:,1]; # number of atoms per each configuration 
    config.a = tm[:,2:(1+dim)]'; # 1st principal vectors per each configuration 
    config.b = tm[:,(2+dim):(1+2*dim)]'; # 2nd principal vectors per each configuration 
    if dim==3
        config.c = tm[:,(2+2*dim):(1+dim*dim)]'; # 2nd principal vectors per each configuration 
    else
        config.c = [];
    end
    if nce>0
        config.e = tm[:,2+dim*dim]; # potential energy per each configuration 
    end

    config.natomall = sum(config.natom);# number of atoms for all configurations
    m  = (1+ncx+ncq+ncv+ncf)*config.natomall;
    n1 = n2+1;
    n2 = n2+m;    
    tm = reshape(tmp[n1:n2], config.natomall, (1+ncx+ncq+ncv+ncf)); 
    config.t = tm[:,1];
    config.x = tm[:,2:(1+ncx)]';
    config.q = tm[:,(2+ncx):(1+ncx+ncq)]';
    config.v = tm[:,(2+ncx+ncq):(1+ncx+ncq+ncv)]';
    config.f = tm[:,(2+ncx+ncq+ncv):(1+ncx+ncq+ncv+ncf)]';

elseif mode==2
    n2 = 6;
    m  = nconfigs;
    n1 = n2+1;
    n2 = n2+m;
    config.natom = reshape(tmp[n1:n2], 1, m); # number of atoms per each configuration  
    config.natomall = sum(config.natom);# number of atoms for all configurations    
    
    # the 1st principal vector of the simulation box
    m  = dim*nconfigs;
    n1 = n2+1;
    n2 = n2+m;
    config.a = reshape(tmp[n1:n2], dim, nconfigs); # the 1st principal vector of the simulation box

    # the 2nd principal vector of the simulation box
    n1 = n2+1;
    n2 = n2+m;
    config.b = reshape(tmp[n1:n2], dim, nconfigs); # the 2nd principal vector of the simulation box

    if dim==3
        # the 3rd principal vector of the simulation box
        n1 = n2+1;
        n2 = n2+m;
        config.c = reshape(tmp[n1:n2], dim, nconfigs); # the 3rd principal vector of the simulation box
    end

    m  = nconfigs;
    n1 = n2+1;
    n2 = n2+m;
    config.e = reshape(tmp[n1:n2], 1, m); # number of atoms per each configuration  
    
    # read atom types for all configurations
    m  = config.natomall;
    n1 = n2+1;
    n2 = n2+m;
    config.t = reshape(int64(tmp[n1:n2]), 1, m);

    # read atom positions for all configurations
    m  = ncx*config.natomall;
    n1 = n2+1;
    n2 = n2+m;
    config.x = reshape(tmp[n1:n2], ncx, config.natomall);

    # read atom charges for all configurations
    m  = ncq*config.natomall;
    n1 = n2+1;
    n2 = n2+m;
    config.q = reshape(tmp[n1:n2], ncq, config.natomall);

    # read atom velocities for all configurations
    m  = ncv*config.natomall;
    n1 = n2+1;
    n2 = n2+m;
    config.v = reshape(tmp[n1:n2], ncv, config.natomall);

    # read atom forces for all configurations
    m  = ncf*config.natomall;
    n1 = n2+1;
    n2 = n2+m;
    config.f = reshape(tmp[n1:n2], ncf, config.natomall);   
end

return config

end
