#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct DomainStruct

    boxlo::Array{Float64,2};
    boxhi::Array{Float64,2};
    boxtilt::Array{Float64,2};
    pbc::Array{Int64,2};
    bcs::Array{Int64,2};
    triclinic::Int64;          

    DomainStruct() = new();
end

function setdomain(boxhi, boxtilt, pbc=[1 1 1], bcs=[1 1 1 1 1 1])

dom = DomainStruct();

reg.boxlo = [0.0 0.0 0.0];
dom.boxhi = boxhi;
dom.boxtilt = boxtilt;
dom.pbc = pbc;
dom.bcs = bcs;

if (dom.boxtilt[1] == 0 && dom.boxtilt[2] == 0 && dom.boxtilt[3] == 0)
    dom.triclinic = 0;
else    
    dom.triclinic = 1;
end

return dom

end

