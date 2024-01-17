#***************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
#  
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#***************************************************************************

function readMDP(filename, mode)

# natom, lattice, pbc, energy, stress, id, type, group, move, atomic number, atomic mass, charge, position, velocity, force

if mode==0 # binary    
    frame = reinterpret(Float64,read(filename));
else # text    
    tmp = Main.DelimitedFiles.readdlm(filename, ',');
    tmp = String.(reduce(vcat, split.(tmp))); # convert to 1D array of strings
    frame = parse.(Float64, tmp);  # convert to 1D array of Float64  
end

config = initializeconfig(0);  

natomall = 0;
nframe = 0;
n2 = 0;
n = length(frame);
while (n2 < n)
    nframe = nframe + 1;
    
    # natom, lattice, pbc, energy, stress, id, type, group, move, atomic number, atomic mass, charge, position, velocity, force
    natom = Int64(frame[n2+1]);  # number of atoms
    ncl = Int64(frame[n2+2]);    # lattice  
    ncp = Int64(frame[n2+3]);    # periodic boundary conditions     
    nce = Int64(frame[n2+4]);    # energy
    ncs = Int64(frame[n2+5]);    # stress 
    nci = Int64(frame[n2+6]);    # atom id
    nct = Int64(frame[n2+7]);    # atom type
    ncg = Int64(frame[n2+8]);    # atom groups
    nco = Int64(frame[n2+9]);    # atom moved/fixed flags
    ncz = Int64(frame[n2+10]);   # atomic numbers
    ncm = Int64(frame[n2+11]);   # atomic masses     
    ncq = Int64(frame[n2+12]);   # atom charges
    ncx = Int64(frame[n2+13]);   # positions
    ncv = Int64(frame[n2+14]);   # atom velocities 
    ncf = Int64(frame[n2+15]);   # atom forces
    n2 = n2+15;
            
    if nframe==1
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
        config.nconfigs = 1; 
        config.dim = ncx;
    end

    config.natom = hcat(config.natom, reshape([natom],1,1))
    
    # lattice, pbc, energy, stress
    n1 = n2 + 1; n2 = n2 + ncl;
    if nframe == 1 & ncl>0
        config.lattice = reshape(frame[n1:n2], (ncl, 1));
    elseif ncl>0
        config.lattice = hcat(config.lattice, reshape(frame[n1:n2], (ncl, 1)));
    end            

    n1 = n2 + 1; n2 = n2 + ncp;
    if nframe == 1 & ncp>0
        config.pbc = reshape(frame[n1:n2], (ncp, 1));
    elseif ncp>0
        config.pbc = hcat(config.pbc, reshape(frame[n1:n2], (ncp, 1)));
    end            

    n1 = n2 + 1; n2 = n2 + nce;
    if nframe == 1 & nce>0
        config.e = reshape(frame[n1:n2], (nce, 1));
    elseif nce>0
        config.e = hcat(config.e, reshape(frame[n1:n2], (nce, 1)));
    end            

    n1 = n2 + 1; n2 = n2 + ncs;
    if nframe == 1 & ncs>0
        config.stress = reshape(frame[n1:n2], (ncs, 1));
    elseif ncs>0
        config.stress = hcat(config.stress, reshape(frame[n1:n2], (ncs, 1)));
    end            
    
    k = nci + nct + ncg + nco + ncz + ncm + ncq + ncx + ncv + ncf;
    n1 = n2 + 1; n2 = n2 + natom*k;    
    tmp = reshape(frame[n1:n2], (k, natom));            

    if nframe == 1 & nci>0
        config.tags = tmp[1:nci,:];
    elseif nci>0
        config.tags = hcat(config.tags, tmp[1:nci,:]);
    end            

    if nframe == 1 & nct>0
        config.t = tmp[(1+nci):(nci+nct),:];               
    elseif nct>0
        config.t = hcat(config.t, tmp[(1+nci):(nci+nct),:]);
    end            

    if nframe == 1 & ncg>0
        config.group = tmp[(1+nci+nct):(nci+nct+ncg),:];
    elseif ncg>0
        config.group = hcat(config.group, tmp[(1+nci+nct):(nci+nct+ncg),:]);
    end            

    if nframe == 1 & nco>0
        config.move = tmp[(1+nci+nct+ncg):(nci+nct+ncg+nco),:];
    elseif nco>0
        config.move = hcat(config.move, tmp[(1+nci+nct+ncg):(nci+nct+ncg+nco),:]);
    end            

    if nframe == 1 & ncz>0
        config.Z = tmp[(1+nci+nct+ncg+nco):(nci+nct+ncg+nco+ncz),:];
    elseif ncz>0
        config.Z = hcat(config.Z, tmp[(1+nci+nct+ncg+nco):(nci+nct+ncg+nco+ncz),:]);
    end            

    if nframe == 1 & ncm>0
        config.mass = tmp[(1+nci+nct+ncg+nco+ncz):(nci+nct+ncg+nco+ncz+ncm),:];
    elseif ncm>0
        config.mass = hcat(config.mass, tmp[(1+nci+nct+ncg+nco+ncz):(nci+nct+ncg+nco+ncz+ncm),:]);
    end            

    if nframe == 1 & ncq>0
        config.q = tmp[(1+nci+nct+ncg+nco+ncz+ncm):(nci+nct+ncg+nco+ncz+ncm+ncq),:];
    elseif ncq>0
        config.q = hcat(config.q, tmp[(1+nci+nct+ncg+nco+ncz+ncm):(nci+nct+ncg+nco+ncz+ncm+ncq),:]);
    end            

    if nframe == 1 & ncx>0
        config.x = tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx),:];
    elseif ncx>0
        config.x = hcat(config.x, tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx),:]);
    end            

    if nframe == 1 & ncv>0
        config.v = tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv),:];
    elseif ncv>0
        config.v = hcat(config.v, tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv),:]);
    end            

    if nframe == 1 & ncf>0
        config.f = tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv+ncf),:];
    elseif ncf>0
        config.f = hcat(config.f, tmp[(1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv+ncf),:]);
    end            
    
    natomall = natomall + natom;
end
config.natomall = natomall;

return config 

end



