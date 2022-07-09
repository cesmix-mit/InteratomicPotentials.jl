#***************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
#  
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#***************************************************************************

function readconfigfile(app)

print("Read configuration from " * app.datapath);
filename = app.datapath * "/" *  app.datafile;
style = lowercase(app.dataformat);

if style == "lammps"     
    mode = app.datamode;
    if mode==0 # binary    
        tmp = reinterpret(Float64,read(filename));
    else # text    
        tmp = Main.DelimitedFiles.readdlm(filename);
    end    
    config = readLAMMPS(tmp);
    config.nconfigs = length(config.natom); # number of configurations
    config.natomall = sum(config.natom);# number of atoms for all configurations    
elseif style == "mdp" 
    mode = app.datamode;
    config = readMDP(filename, mode);
    dim = config.dim;    
    if size(config.lattice,1) == dim*dim    
        config.a = config.lattice[1:dim,:];
        config.b = config.lattice[(dim+1):2*dim,:];
        if dim==3
            config.c = config.lattice[(2*dim+1):3*dim,:];
        end
    end
elseif style == "extxyz" 
    config = readEXTXYZ(filename, app.atomspecies);    
    dim = config.dim;    
    if size(config.lattice,1) == dim*dim    
        config.a = config.lattice[1:dim,:];
        config.b = config.lattice[(dim+1):2*dim,:];
        if dim==3
            config.c = config.lattice[(2*dim+1):3*dim,:];
        end
    end
elseif style == "json" 
    config = readJSON(filename, app.atomspecies);    
    dim = config.dim;    
    if size(config.lattice,1) == dim*dim    
        config.a = config.lattice[1:dim,:];
        config.b = config.lattice[(dim+1):2*dim,:];
        if dim==3
            config.c = config.lattice[(2*dim+1):3*dim,:];
        end
    end
end

return config

end

