#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function domainmaps(a::Array{Float64,1}, b::Array{Float64,1})
    # a, b, c are three vectors describing the simulation box
    # a = [a1, a2, a3]
    # b = [b1, b2, b3]
    # c = [c1, c2, c3]
    
    R2B = [a b]; # transform the unit square to the simulation box (a, b)  
    B2R = inv(R2B); # transform the simulation box (a, b) to the unit square 
    
    return B2R, R2B
    
end

function domainmaps(a::Array{Float64,1}, b::Array{Float64,1}, c::Array{Float64,1})
    # a, b, c are three vectors describing the simulation box
    # a = [a1, a2, a3]
    # b = [b1, b2, b3]
    # c = [c1, c2, c3]
    
    R2B = [a b c]; # transform the unit cube to the simulation box (a, b, c)  
    B2R = inv(R2B); # transform the simulation box (a, b, c) to the unit cube
    
    return B2R, R2B
    
end
    
