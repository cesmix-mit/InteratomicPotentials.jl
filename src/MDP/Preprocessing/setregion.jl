#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct RegStruct

    boxlo::Array{Float64,2};
    boxhi::Array{Float64,2};
    boxtilt::Array{Float64,2};
    triclinic::Int64;          

    RegStruct() = new();
end

function setregion(boxhi, boxtilt=[0.0 0.0 0.0])

reg = RegStruct();

reg.boxlo = [0.0 0.0 0.0];
reg.boxhi = boxhi;
reg.boxtilt = boxtilt;

if (reg.boxtilt[1] == 0 && reg.boxtilt[2] == 0 && reg.boxtilt[3] == 0)
    reg.triclinic = 0;
else    
    reg.triclinic = 1;
end

return reg

end

