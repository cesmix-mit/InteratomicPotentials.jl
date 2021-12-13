#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct LatticeStruct

    origin::Array{Float64,2};
    orientx::Array{Float64,2};
    orienty::Array{Float64,2};
    orientz::Array{Float64,2};
    spacing::Array{Float64,2};
    a1::Array{Float64,2};
    a2::Array{Float64,2};
    a3::Array{Float64,2};
    basis::Array{Float64,2};
    type::Array{Int64,2};    
    style::Int64;         
    spaceflag::Int64;         
    scale::Float64; 
    LatticeStruct() = new();
end

function setlattice(stylein, scale, orientx=[1.0 0.0 0.0], orienty=[0.0 1.0 0.0], orientz=[0.0 0.0 1.0], spacing=[1.0 1.0 1.0], type=nothing, a1=[1.0 0.0 0.0], a2=[0.0 1.0 0.0], a3=[0.0 0.0 1.0], basis=nothing)

lat = LatticeStruct();

lat.scale = scale;
lat.spaceflag = 0;

if (lowercase(stylein) == "none") 
    lat.style = 0;
elseif (lowercase(stylein) == "sc") 
    lat.style = 1;    
elseif (lowercase(stylein) == "bcc") 
    lat.style = 2;
elseif (lowercase(stylein) == "fcc") 
    lat.style = 3;
elseif (lowercase(stylein) == "hcp") 
    lat.style = 4;
elseif (lowercase(stylein) == "diamond") 
    lat.style = 5;
elseif (lowercase(stylein) == "sq") 
    lat.style = 6;
elseif (lowercase(stylein) == "sq2") 
    lat.style = 7;
elseif (lowercase(stylein) == "hex") 
    lat.style = 8;        
elseif (lowercase(stylein) == "custom") 
    lat.style = 9;        
else
    error("Invalid lattice style");
end

lat.origin = [0.0 0.0 0.0];
lat.orientx = orientx;
lat.orienty = orienty;
lat.orientz = orientz;
lat.spacing = spacing;

if (lowercase(stylein) == "custom") 
    lat.a1 = a1;
    lat.a2 = a2;
    lat.a3 = a3;
    lat.basis = basis;    
else
    lat.a1 = [1.0 0.0 0.0];
    lat.a2 = [0.0 1.0 0.0];
    lat.a3 = [0.0 0.0 1.0];
    
    if (lowercase(stylein) === "hex") 
        lat.a2[2] = sqrt(3.0);
    end
    if (lowercase(stylein) === "hcp") 
        lat.a2[2] = sqrt(3.0);
        lat.a3[3] = sqrt(8.0/3.0);
    end
    
    if (lowercase(stylein) === "sc") 
        lat.basis = [0 0 0];        
    elseif (lowercase(stylein) === "bcc") 
        lat.basis = [0 0 0; 0.5 0.5 0.5];
    elseif (lowercase(stylein) === "fcc") 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5];
    elseif (lowercase(stylein) === "hcp") 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 5.0/6.0 0.5; 0.0 1.0/3.0 0.5];
    elseif (lowercase(stylein) === "diamond") 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5; 0.25 0.25 0.25; 0.25 0.75 0.75; 0.75 0.25 0.75; 0.75 0.75 0.25];        
    elseif (lowercase(stylein) === "sq") 
        lat.basis = [0 0 0];        
    elseif (lowercase(stylein) === "sq2") 
        lat.basis = [0 0 0; 0.5 0.5 0.0];
    elseif (lowercase(stylein) === "hex") 
        lat.basis = [0 0 0; 0.5 0.5 0.0];    
    end
end

nbasis = size(lat.basis)[1];

if type===nothing
    lat.type = Int64.(ones(1,nbasis));
else
    lat.type = type;
end

return lat

end

