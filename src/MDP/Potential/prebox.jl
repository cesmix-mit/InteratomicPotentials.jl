#****************************************************************************                               
#                     Molecular Dynamics Potentials [MDP]
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen [cuongng@mit.edu, exapde@gmail.com]
#****************************************************************************

function prebox(pbc::Array{Int64,1}, rmax::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1}, c=nothing)

if c===nothing
    B2R, R2B = domainmaps(a, b);

    v,w = boundingbox(pbc, rmax, a, b);
    vr = B2R*v; wr = B2R*w;

    pimages = periodicimages(pbc, a, b);    
else
    B2R, R2B = domainmaps(a, b, c);

    v,w = boundingbox(pbc, rmax, a, b, c);
    vr = B2R*v; wr = B2R*w;

    pimages = periodicimages(pbc, a, b, c);
end

return B2R, R2B, v, w, vr, wr, pimages

end

