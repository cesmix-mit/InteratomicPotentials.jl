function parsecomment(s)

p = split(s, '=', keepempty=false);
n = length(p);

keys = Array{Any}(undef,n-1);
values = Array{Any}(undef,n-1);
for i = 1:n-1
    ii = findlast("\"", p[i]);
    if ii === nothing
        q = split(p[i], ' ', keepempty=false);
        if length(q) == 1
            keys[i] = q[1];
        elseif length(q)==2
            keys[i] = q[2];
            values[i-1] = q[1];
        end
    else
        t = split(p[i], '"', keepempty=false);
        q = split(t[1], ' ', keepempty=false);
        values[i-1] = map(x -> parse(Float64, x), q);
        keys[i] = t[2][2:end];    
    end
end
t = split(p[n], '"', keepempty=false);
values[n-1] = t[1];

# handle energy
for i = 1:n-1    
    if lowercase(string(keys[i])) == "energy"                
        values[i] = parse(Float64, values[i]);        
    end
end

# handle pbc
for i = 1:n-1    
    if lowercase(string(keys[i])) == "pbc"        
        values[i] = replace(values[i], "T" => "1");
        values[i] = replace(values[i], "F" => "0");        
        q = split(values[i], ' ', keepempty=false);
        values[i] = map(x -> parse(Int64, x), q);        
    end
end

# handle properties
indp = 0;
for i = 1:n-1    
    if lowercase(string(keys[i])) == "properties"        
        values[i] = split(values[i], ':', keepempty=false);
        indp = i;
        break;        
    end
end

return keys, values, indp

end

function readEXTXYZ(filename, species)

config = initializeconfig(0);  

f = readlines(filename);
n = length(f);

natomall = 0;
nframe = 0;
m = 0;
while (m<n)
    m = m + 1;
    nframe = nframe + 1          
    natom = parse(Int64, f[m])    
    config.natom = hcat(config.natom, reshape([natom],1,1))

    m = m + 1
    keys, values, indp = parsecomment(f[m]);    
    for i = 1:length(keys)
        mystr = lowercase(string(keys[i]));
        if mystr == "lattice"
            d = reshape(values[i], length(values[i]), 1);
            if nframe == 1
                config.lattice = d;
            else
                config.lattice = hcat(config.lattice, d);
            end            
            config.ncl = length(values[i]);
        elseif mystr == "energy"            
            d = reshape([values[i]], 1, 1);
            if nframe == 1
                config.e = d;
            else    
                config.e = hcat(config.e, d);
            end            
            config.nce = 1;
        elseif mystr == "stress"
            d = reshape(values[i], length(values[i]), 1);
            if nframe == 1
                config.stress = d;
            else    
                config.stress = hcat(config.stress, d);
            end            
            config.ncs = length(values[i]);
        elseif mystr == "pbc"
            d = reshape(values[i], length(values[i]), 1);
            if nframe == 1
                config.pbc = d;
            else    
                config.pbc = hcat(config.pbc, d);
            end            
            config.ncp = length(values[i]);
        end
    end
        
    s = String.(reduce(vcat, split.(f[(m+1):(m+natom)])))
    s = reshape(s, (Int64(length(s)/natom), natom))
    m = m + natom    
    
    # display(keys)
    # display(values)
    # display(indp)
    # display(values[indp])
    # handle properties
    prop = reshape(values[indp],(3, Int64(length(values[indp])/3)));
    nfield = size(prop,2);            
    fieldtype = String.(prop[2,:]);
    fieldlength = parse.(Int64, prop[3,:]);
    fieldlength = Int64.(cumsum(fieldlength[:]));    
    
    frame = Array{Any}(undef,nfield);
    for j = 1:nfield        
        if j == 1
            k1 = 1;
            k2 = fieldlength[1];
        else
            k1 = fieldlength[j-1]+1;
            k2 = fieldlength[j];
        end
        if lowercase(fieldtype[j]) == "s" || lowercase(fieldtype[j]) == "l"               
            frame[j] = s[k1,:]; 
            b = ones(natom,1);
            if lowercase(fieldtype[j]) == "l"                              
                ind = findall(frame[j][:] .== "F");
                b[ind] .= 0;            
                frame[j] = b;
            else
                for k = 1:length(species)
                    ind = findall(frame[j][:] .== species[k]);
                    b[ind] .= k;
                end              
                frame[j] = b;                                  
            end           
        else
           frame[j] = parse.(Float64, s[k1:k2,:]); 
        end        
    end    

    for j = 1:nfield
        mystr = lowercase(prop[1,j]);
        if mystr == "species"
            if nframe == 1
                config.t = frame[j];
            else    
                config.t = hcat(config.t, frame[j]);
            end
            config.nct = size(config.t,1);                    
        elseif mystr == "pos"
            if nframe == 1
                config.x = frame[j];
            else    
                config.x = hcat(config.x, frame[j]);
            end
            config.ncx = size(config.x,1);
            config.dim = size(config.x,1);
        elseif mystr == "vel"
            if nframe == 1
                config.v = frame[j];
            else    
                config.v = hcat(config.v, frame[j]);
            end
            config.ncv = size(config.v,1);
        elseif mystr == "forces"
            if nframe == 1
                config.f = frame[j];
            else    
                config.f = hcat(config.f, frame[j]);
            end
            config.ncf = size(config.f,1);            
        elseif mystr == "charge"
            if nframe == 1
                config.q = frame[j];
            else    
                config.q = hcat(config.q, frame[j]);
            end
            config.ncq = size(config.q,1);        
        elseif mystr == "mass"
            if nframe == 1
                config.mass = frame[j];
            else    
                config.mass = hcat(config.mass, frame[j]);
            end
            config.ncm = size(config.mass,1);                                 
        elseif mystr == "tags"
            if nframe == 1
                config.tags = frame[j];
            else    
                config.tags = hcat(config.tags, frame[j]);
            end
            config.nci = size(config.tags,1);                    
        elseif mystr == "groups"
            if nframe == 1
                config.group = frame[j];
            else    
                config.group = hcat(config.group, frame[j]);
            end
            config.ncg = size(config.group,1);                                                                                                                                                              
        elseif mystr == "move_mask"
            if nframe == 1
                config.move = frame[j];
            else    
                config.move = hcat(config.move, frame[j]);
            end
            config.nco = size(config.move,1);                                                                                                                                                             
        elseif mystr == "z"
            if nframe == 1
                config.Z = frame[j];
            else    
                config.Z = hcat(config.Z, frame[j]);
            end
            config.ncz = size(config.Z,1);                                                                                                                                                                          
        end
    end
end

return config;

end

