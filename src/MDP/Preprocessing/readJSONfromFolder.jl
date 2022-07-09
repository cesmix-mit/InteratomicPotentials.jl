#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function readJSONfromFolder(datapath, species)

    config = initializeconfig(0);  

    filenames = readdir(datapath);
    n = length(filenames);
    nframe = 0;
    for i = 1:n
        fname = datapath * "/" * filenames[i];
        mystr = read(fname, String);
        i1 = findfirst("{", mystr);
        mystr = mystr[i1[1]:end];
        dict2 = JSON.parse(mystr);
        Dataset = dict2["Dataset"]["Data"];

        #println(Dataset)
        for j = 1:length(Dataset)
            Data = Dataset[j];

            nframe = nframe + 1;

            if haskey(Data, "NumAtoms")
                config.natom = hcat(config.natom, reshape([Data["NumAtoms"]],1,1))
                natom = Data["NumAtoms"];       
            end    
            if haskey(Data, "Lattice")
                d = hcat(Data["Lattice"]...);
                d = reshape(d, length(d), 1);
                if nframe == 1
                    config.lattice = d;
                else                        
                    config.lattice = hcat(config.lattice, d);
                end
                config.ncl = size(config.lattice,1);
            end
            if haskey(Data, "Stress")
                d = hcat(Data["Stress"]...);
                d = reshape(d, length(d), 1);
                if nframe == 1
                    config.stress = d;
                else    
                    config.stress = hcat(config.stress, d);
                end
                config.ncs = size(config.stress,1);
            end
            if haskey(Data, "Energy")
                d = hcat(Data["Energy"]...);
                d = reshape(d, length(d), 1);
                if nframe == 1
                    config.e = d;
                else    
                    config.e = hcat(config.e, d);
                end
                config.nce = size(config.e,1);
            end
            if haskey(Data, "Positions")
                d = hcat(Data["Positions"]...);
                if nframe == 1
                    config.x = d;
                else    
                    config.x = hcat(config.x, d);
                end
                config.ncx = size(config.x,1);
                config.dim = size(config.x,1);
            end
            if haskey(Data, "Forces")
                d = hcat(Data["Forces"]...);
                if nframe == 1
                    config.f = d;
                else    
                    config.f = hcat(config.f, d);
                end
                config.ncf = size(config.f,1);
            end
            if haskey(Data, "Velocities")
                d = hcat(Data["Velocities"]...);
                if nframe == 1
                    config.v = d;
                else    
                    config.v = hcat(config.v, d);
                end
                config.ncf = size(config.v,1);
            end
            if haskey(Data, "Charges")
                d = hcat(Data["Charges"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.q = d;
                else    
                    config.q = hcat(config.q, d);
                end
                config.ncq = size(config.q,1);
            end
            if haskey(Data, "Masses")
                d = hcat(Data["Masses"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.mass = d;
                else    
                    config.mass = hcat(config.mass, d);
                end
                config.ncm = size(config.mass,1);
            end
            if haskey(Data, "Groups")
                d = hcat(Data["Groups"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.group = d;
                else    
                    config.group = hcat(config.group, d);
                end
                config.ncg = size(config.group,1);
            end
            if haskey(Data, "Moves")
                d = hcat(Data["Moves"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.move = d;
                else    
                    config.move = hcat(config.move, d);
                end
                config.nco = size(config.move,1);
            end
            if haskey(Data, "Tags")
                d = hcat(Data["Tags"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.tags = d;
                else    
                    config.tags = hcat(config.tags, d);
                end
                config.nct = size(config.tags,1);
            end
            if haskey(Data, "Z")
                d = hcat(Data["Z"]...);
                d = reshape(d, 1, length(d));
                if nframe == 1
                    config.Z = d;
                else    
                    config.Z = hcat(config.Z, d);
                end
                config.ncz = size(config.Z,1);
            end
            if haskey(Data, "AtomTypes")
                d = hcat(Data["AtomTypes"]...);
                d = reshape(d, 1, length(d));
                b = ones(1, length(d));
                for k=1:length(species)
                    ind = findall(d[:] .== species[k]);
                    b[ind] .= k;
                end                                
                if nframe == 1
                    config.t = b;
                else    
                    config.t = hcat(config.t, b);
                end
                config.nct = size(config.t,1);
            end
        end
    end
    config.nconfigs = nframe;
    config.natomall = sum(config.natom); 

    return config

end
